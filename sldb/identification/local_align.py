import itertools
import multiprocessing as mp
import re
import subprocess
import shlex
from sqlalchemy import func

import dnautils
import sldb.util.funcs as funcs
from sldb.identification import AlignmentException
from sldb.identification.j_genes import JGermlines
from sldb.identification.v_genes import VGermlines
from sldb.identification.vdj_sequence import VDJSequence
from sldb.common.models import (CDR3_OFFSET, DuplicateSequence, NoResult,
                                Sample, Sequence, serialize_gaps)
import sldb.util.lookups as lookups
import sldb.util.concurrent as concurrent

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
                 min_similarity, max_deletions, max_insertions):
        self.complete_queue = complete_queue
        self.v_germlines = v_germlines
        self.j_germlines = j_germlines
        self.align_path = align_path
        self.min_similarity = min_similarity
        self.max_deletions = max_deletions
        self.max_insertions = max_insertions

        self.first_alleles = {
            name: v.sequence for name, v in self.v_germlines.iteritems()
            if int(name.split('*', 1)[1]) == 1
        }

    def do_task(self, args):
        record = {
            'seq_id': args['seq_ids'][0],
            'sample_id': args['sample_id'],
            'subject_id': args['subject_id'],
        }
        # Find best aligned first allele
        v_align = self.align_seq_to_germs(args['seq'], self.first_alleles)

        if v_align is None or not self.alignment_passes(v_align['germ'],
                                                        v_align['seq']):
            return

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

        pre_cdr3_germ = v_align['germ'][:cdr3_start].replace(
            GAP_PLACEHOLDER, '')
        pre_cdr3_seq = v_align['seq'][:cdr3_start].replace(
            GAP_PLACEHOLDER, '')
        v_length = len(pre_cdr3_seq)
        v_match = v_length - dnautils.hamming(pre_cdr3_germ, pre_cdr3_seq)
        if v_match / float(v_length) < self.min_similarity:
            return

        # NOTE: This doesn't look for a streak like VDJSequence
        record.update({
            'v_match': v_match,
            'v_length': v_length,
            'v_mutation_fraction': v_match / float(v_length),
            'pre_cdr3_match': v_match,
            'pre_cdr3_length': v_length,
            'probable_indel_or_misalign': False,
            'insertions': gap_positions(v_align['germ'][:cdr3_start]),
            'deletions': gap_positions(v_align['seq'][:cdr3_start]),
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
            return

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
            return

        final_germ = ''.join((
            final_germ[:cdr3_start],
            '-' * (cdr3_end - cdr3_start),
            final_germ[cdr3_end:]
        ))

        stop = lookups.has_stop(final_seq)
        in_frame = cdr3_start % 3 == 0 and cdr3_end % 3 == 0

        post_cdr3_germ = final_germ[-self.j_germlines.upstream_of_cdr3:]
        post_cdr3_seq = final_seq[-self.j_germlines.upstream_of_cdr3:]
        post_cdr3_length = len(post_cdr3_seq)
        j_match = post_cdr3_length - dnautils.hamming(
            post_cdr3_germ,
            post_cdr3_seq
        )

        record.update({
            'v_gene': v_align['germ_name'],
            'j_gene': j_align['germ_name'],

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

            # NOTE: These  may not be correct
            'paired': True,

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
            '{} --match 2 --mismatch -2 --gapopen -5 --gapextend -2 '
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

        if task['type'] == 'Sequence':
            seq = session.query(Sequence).filter(
                Sequence.sample_id == task['record']['sample_id'],
                Sequence.seq_id == task['record']['seq_id']
            ).one()
            for key, value in task['record'].iteritems():
                setattr(seq, key, value)
        else:
            try:
                new_seq = Sequence(**task['record'])
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
            except ValueError as e:
                pass

        session.commit()


def run_fix_sequences(session, args):
    v_germlines = VGermlines(args.v_germlines)
    j_germlines = JGermlines(args.j_germlines, args.upstream_of_cdr3, 0, 0)

    mutation_cache = {}
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
        print ('Creating task queue for sample {}; '
               '{} indels, {} noresults').format(sample_id, indels.count(),
                                                 noresults.count())

        # Get information for the V-ties
        avg_mut, avg_len = session.query(
            Sample.v_ties_mutations,
            Sample.v_ties_len
        ).filter(Sequence.sample_id == sample_id).first()

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
                args.max_insertions)
            )

        tasks.start(block=False)

        process_completes(session, complete_queue, workers)
