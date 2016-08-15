import itertools
import multiprocessing as mp
import re
import subprocess
import shlex
from sqlalchemy import desc

from Bio.Seq import Seq
import dnautils
from immunedb.identification.j_genes import JGermlines
from immunedb.identification.v_genes import VGermlines
from immunedb.common.models import (CDR3_OFFSET, DuplicateSequence, NoResult,
                                    Sample, Sequence)
import immunedb.util.concurrent as concurrent
import immunedb.util.funcs as funcs
from immunedb.util.log import logger
import immunedb.util.lookups as lookups

GAP_PLACEHOLDER = '^'


def gaps_before(gaps, pos):
    return sum((e[1] for e in gaps if e[0] < pos))


def gap_positions(seq):
    gaps = []
    for diff in re.finditer('[-]+', seq):
        start, end = diff.span()
        gaps.append((start, end - start))
    return gaps


class LocalAlignmentWorker(concurrent.Worker):
    def __init__(self, complete_queue, v_germlines, j_germlines, align_path,
                 min_similarity, max_deletions, max_insertions, max_padding):
        self.complete_queue = complete_queue
        self.v_germlines = v_germlines
        self.j_germlines = j_germlines
        self.align_path = align_path
        self.min_similarity = min_similarity
        self.max_deletions = max_deletions
        self.max_insertions = max_insertions
        self.max_padding = max_padding

        self.first_alleles = {
            name: v.sequence
            for name, v in self.v_germlines.alignments.iteritems()
            if int(name.split('*', 1)[1]) == 1
        }

    def do_task(self, args, rc=False):
        if rc:
            args['seq'] = str(Seq(args['seq']).reverse_complement())

        record = {
            'seq_id': args['seq_ids'][0],
            'sample_id': args['sample_id'],
            'subject_id': args['subject_id'],
        }
        # Find best aligned first allele
        v_align = self.align_seq_to_germs(args['seq'], self.first_alleles)

        if (v_align is None or
                not self.alignment_passes(v_align['germ'], v_align['seq']) or
                (self.max_padding is not None and
                    v_align['seq_offset'] > self.max_padding)):
            return None if rc else self.do_task(args, True)

        v_name = v_align['germ_name'].split('*', 1)[0]
        v_ties = {
            '|'.join(name): v for name, v in self.v_germlines.all_ties(
                args['avg_len'], args['avg_mut'], cutoff=False
            ).iteritems() if v_name in '|'.join(name)
        }

        v_align = self.align_seq_to_germs(
            v_align['seq'].replace('-', ''), v_ties
        )

        germ_insertions = gap_positions(v_align['germ'])
        imgt_gaps = [
            i + gaps_before(germ_insertions, i)
            for i, c in enumerate(v_ties[v_align['germ_name']]) if c == '-'
        ]
        for gap in imgt_gaps:
            v_align['germ'] = ''.join((
                v_align['germ'][:gap],
                GAP_PLACEHOLDER,
                v_align['germ'][gap:]
            ))
            v_align['seq'] = ''.join((
                v_align['seq'][:gap],
                GAP_PLACEHOLDER,
                v_align['seq'][gap:]
            ))
        cdr3_start = 0
        count = 0
        for i, c in enumerate(v_align['germ']):
            if c != '-':
                count += 1
                if count > CDR3_OFFSET:
                    cdr3_start = i
                    break
        v_align['seq'] = ''.join((
            v_align['seq'][:cdr3_start],
            v_align['seq'][cdr3_start:].replace('-', '')
        ))

        pre_cdr3_germ = v_align['germ'][:cdr3_start].replace(
            GAP_PLACEHOLDER, '')
        pre_cdr3_seq = v_align['seq'][:cdr3_start].replace(
            GAP_PLACEHOLDER, '')
        v_length = len(pre_cdr3_seq)
        v_match = v_length - dnautils.hamming(pre_cdr3_germ, pre_cdr3_seq)
        if v_match / float(v_length) < self.min_similarity:
            return None if rc else self.do_task(args, True)

        # NOTE: This doesn't look for a streak like VDJSequence
        record.update({
            'locally_aligned': True,
            'v_match': v_match,
            'v_length': v_length,
            'v_mutation_fraction': v_match / float(v_length),
            'pre_cdr3_match': v_match,
            'pre_cdr3_length': v_length,
            'probable_indel_or_misalign': False,
            'insertions': set(gap_positions(v_align['germ'][:cdr3_start])),
            'deletions': set(gap_positions(v_align['seq'][:cdr3_start])),
        })

        v_align['germ'] = v_align['germ'].replace(GAP_PLACEHOLDER, '-')
        v_align['seq'] = v_align['seq'].replace(GAP_PLACEHOLDER, '-')
        record.update({
            'num_gaps': v_align['seq'][:cdr3_start].count('-'),
            'pad_length': v_align['seq_offset']
        })

        j_align = self.align_seq_to_germs(
            v_align['seq'][cdr3_start:], self.j_germlines,
        )

        if j_align is None:
            return None if rc else self.do_task(args, True)

        final_germ = ''.join((
            v_align['germ'][:-len(j_align['germ'])],
            j_align['germ']
        ))
        final_seq = ''.join((
            v_align['seq'][:-len(j_align['seq'])],
            j_align['seq']
        ))
        final_germ = final_germ.rstrip('-')
        final_seq = final_seq[:len(final_germ)]
        cdr3_end = len(final_seq) - self.j_germlines.upstream_of_cdr3

        if cdr3_end - cdr3_start < 3:
            return None if rc else self.do_task(args, True)

        final_germ = ''.join((
            final_germ[:cdr3_start],
            '-' * (cdr3_end - cdr3_start),
            final_germ[cdr3_end:]
        ))
        record['insertions'].update([
            (p[0] + cdr3_end, p[1])
            for p in gap_positions(final_germ[cdr3_end:])
        ])
        record['deletions'].update([(
            p[0] + cdr3_end, p[1])
            for p in gap_positions(final_seq[cdr3_end:])
        ])

        stop = lookups.has_stop(final_seq)
        in_frame = cdr3_start % 3 == 0 and cdr3_end % 3 == 0

        post_cdr3_germ = final_germ[-self.j_germlines.upstream_of_cdr3:]
        post_cdr3_seq = final_seq[-self.j_germlines.upstream_of_cdr3:]
        post_cdr3_length = len(post_cdr3_seq)
        j_match = post_cdr3_length - dnautils.hamming(
            post_cdr3_germ,
            post_cdr3_seq
        )

        v_prefix = v_align['germ_name'][:4]
        j_prefix = j_align['germ_name'][:4]
        record.update({
            'v_gene': funcs.format_ties(
                v_align['germ_name'].split('|'), v_prefix, strip_alleles=True),
            'j_gene': funcs.format_ties(
                j_align['germ_name'].split('|'), j_prefix, strip_alleles=True),

            'cdr3_num_nts': cdr3_end - cdr3_start,
            'cdr3_nt': final_seq[cdr3_start:cdr3_end],
            'cdr3_aa': lookups.aas_from_nts(final_seq[cdr3_start:cdr3_end]),

            'post_cdr3_match': j_match,
            'post_cdr3_length': post_cdr3_length,
            'j_match': j_match,
            'j_length': post_cdr3_length,

            'sequence': final_seq,
            'partial': record['pad_length'] > 0,
            # TODO: Quality
            'germline': final_germ,

            'stop': stop,
            'in_frame': in_frame,
            'functional': stop and in_frame,
            'copy_number': len(args['seq_ids'])
        })

        self.complete_queue.put({
            'type': args['type'],
            'duplicates': args['seq_ids'][1:],
            'record': record
        })

    def cleanup(self):
        self.complete_queue.put(None)

    def align_seq_to_germs(self, seq, germs):
        stdin = []
        for g_name, g_seq in germs.iteritems():
            stdin.append('>{}\n{}\n'.format(g_name, g_seq.replace('-', '')))
            stdin.append('>query\n{}\n'.format(seq.lstrip('N')))

        cmd = (
            '{} --match 2 --mismatch -2 --gapopen -10 --gapextend -5 '
            '--wildcard N 2 --printfasta --printscores --freestartgap '
            '--freeendgap --file -'
        ).format(self.align_path)
        proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, stdin=subprocess.PIPE)
        output, err = proc.communicate(''.join(stdin))
        regex = (
            r'(?P<germ_name>.+)\n'
            r'(?P<germ_padding>[-]*)(?P<germ>[ATCGN-]+)\n'
            r'.+\n'
            r'(?P<seq_padding>[-]*)(?P<seq>[ATCGN-]+)\n'
            r'score: (?P<score>\d+)\n'
        )

        best = None
        for match in re.finditer(regex, output, re.MULTILINE):
            match = match.groupdict()
            match['germ_offset'] = len(match['germ_padding'])
            match['seq_offset'] = len(match['seq_padding'])
            match['score'] = int(match['score'])
            if best is None or match['score'] > best['score']:
                best = match

        if best is None:
            return None
        if best['germ_offset'] > 0:
            best['seq'] = best['seq'][best['germ_offset']:]
        elif best['seq_offset'] > 0:
            best['seq'] = ('N' * best['seq_offset']) + best['seq']

        return best

    def alignment_passes(self, germ, seq):
        if len(gap_positions(seq.strip('-'))) > self.max_deletions:
            return False
        if len(gap_positions(germ.strip('-'))) > self.max_insertions:
            return False
        return True


def process_completes(session, complete_queue, num_workers):
    stops = 0
    while True:
        task = complete_queue.get()
        if task is None:
            stops += 1
            if stops >= num_workers:
                break
            continue

        try:
            if task['type'] == 'Sequence':
                seq = session.query(Sequence).filter(
                    Sequence.sample_id == task['record']['sample_id'],
                    Sequence.seq_id == task['record']['seq_id']
                ).one()
                for key, value in task['record'].iteritems():
                    setattr(seq, key, value)
            else:
                new_seq = Sequence(**task['record'])
                new_seq.locally_aligned = True
                session.add(new_seq)
                # Delete primary NoResult
                session.query(NoResult).filter(
                    NoResult.sample_id == task['record']['sample_id'],
                    NoResult.seq_id == task['record']['seq_id']
                ).delete()

                # Move duplicates to DuplicateSequence
                if len(task['duplicates']) > 0:
                    session.flush()
                    nores = session.query(NoResult).filter(
                        NoResult.sample_id == task['record']['sample_id'],
                        NoResult.seq_id.in_(task['duplicates'])
                    )
                    for old_nores in nores:
                        session.add(DuplicateSequence(
                            seq_id=old_nores.seq_id,
                            sample_id=old_nores.sample_id,
                            duplicate_seq=new_seq))
                        session.delete(old_nores)
        except ValueError:
            pass

        session.commit()


def remove_duplicates(session, sample_id):
    seqs = session.query(
        Sequence.ai, Sequence.seq_id, Sequence.v_gene, Sequence.j_gene,
        Sequence.cdr3_num_nts, Sequence.copy_number, Sequence.sequence
    ).filter(
        Sequence.locally_aligned.is_(True),
        Sequence.sample_id == sample_id
    )

    for seq in seqs:
        potential_collapse = session.query(
            Sequence.ai, Sequence.sequence
        ).filter(
            Sequence.sample_id == sample_id,
            Sequence.v_gene == seq.v_gene,
            Sequence.j_gene == seq.j_gene,
            Sequence.cdr3_num_nts == seq.cdr3_num_nts,
        ).order_by(desc(Sequence.copy_number))

        for other_seq in potential_collapse:
            if (other_seq.ai == seq.ai or
                    len(other_seq.sequence) != len(seq.sequence)):
                continue

            if dnautils.equal(other_seq.sequence, seq.sequence):
                session.query(Sequence).filter(
                    Sequence.ai == other_seq.ai
                ).one().copy_number += seq.copy_number
                session.query(DuplicateSequence).filter(
                    DuplicateSequence.duplicate_seq_ai == seq.ai
                ).update({
                    'duplicate_seq_ai': other_seq.ai,
                }, synchronize_session=False)
                session.add(DuplicateSequence(
                    seq_id=seq.seq_id,
                    duplicate_seq_ai=other_seq.ai,
                    sample_id=sample_id
                ))
                session.query(Sequence).filter(Sequence.ai == seq.ai).delete()
                break
    session.commit()


def run_fix_sequences(session, args):
    v_germlines = VGermlines(args.v_germlines)
    j_germlines = JGermlines(args.j_germlines, args.upstream_of_cdr3, 0, 0)

    for (sample_id, subject_id) in session.query(
            Sample.id, Sample.subject_id).order_by(Sample.id):
        # Get all the indels that were identified (poorly)
        indels = session.query(Sequence).filter(
            Sequence.sample_id == sample_id,
            Sequence.probable_indel_or_misalign == 1
        )
        # Get the sequences that were not identifiable
        noresults = session.query(NoResult).filter(
            NoResult.sample_id == sample_id
        )
        logger.info('Creating task queue for sample {}; '
                    '{} indels, {} noresults'.format(
                        sample_id, indels.count(), noresults.count()))

        # Get information for the V-ties
        avg_mut, avg_len = session.query(
            Sample.v_ties_mutations,
            Sample.v_ties_len
        ).filter(Sample.id == sample_id).one()

        # Get all the unique sequences
        uniques = {}
        for seq in itertools.chain(indels, noresults):
            session.expunge(seq)
            if seq.sequence not in uniques:
                uniques[seq.sequence] = {
                    'type': type(seq).__name__,
                    'sample_id': seq.sample_id,
                    'subject_id': subject_id,
                    'seq_ids': [],
                    'seq': seq.sequence.replace('-', '').strip('N'),
                    'avg_mut': avg_mut,
                    'avg_len': avg_len
                }
            uniques[seq.sequence]['seq_ids'].append(seq.seq_id)

        tasks = concurrent.TaskQueue()
        tasks.add_tasks(uniques.values())

        workers = min(args.nproc, tasks.num_tasks)
        complete_queue = mp.Queue()

        for i in range(0, workers):
            tasks.add_worker(LocalAlignmentWorker(
                complete_queue,
                v_germlines,
                j_germlines,
                args.align_path,
                args.min_similarity / 100.0,
                args.max_deletions,
                args.max_insertions,
                args.max_padding)
            )

        tasks.start(block=False)

        process_completes(session, complete_queue, workers)
        remove_duplicates(session, sample_id)
