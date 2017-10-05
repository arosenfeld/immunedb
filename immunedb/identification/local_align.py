import csv
import os
import re
import subprocess
import shlex
import StringIO

from sqlalchemy import desc

import dnautils

import immunedb.common.config as config
from immunedb.identification import add_as_sequence
from immunedb.identification.vdj_sequence import VDJAlignment, VDJSequence
from immunedb.identification.genes import GeneName, JGermlines, VGermlines
from immunedb.identification.identify import IdentificationProps
from immunedb.common.models import (CDR3_OFFSET, DuplicateSequence, NoResult,
                                    Sample, Sequence, serialize_gaps)
import immunedb.util.concurrent as concurrent
from immunedb.util.funcs import format_ties, periodic_commit
import immunedb.util.lookups as lookups
from immunedb.util.log import logger


GAP_PLACEHOLDER = '.'


def get_fasta(sequences):
    return '\n'.join(['>{}\n{}'.format(k, v.replace('-', ''))
                     for k, v in sequences.iteritems()]) + '\n'


def gaps_before(gaps, pos):
    return sum((e[1] for e in gaps if e[0] < pos))


def gap_positions(seq, char='-'):
    gaps = []
    for diff in re.finditer('[{}]+'.format(char), seq):
        start, end = diff.span()
        gaps.append((start, end - start))
    return gaps


def build_index(germlines, path):
    sequences = get_fasta(germlines)
    cmd = ('bowtie2-build -c "{}" {}').format(sequences, path)
    proc = subprocess.Popen(shlex.split(cmd), stderr=subprocess.PIPE,
                            stdout=subprocess.PIPE)

    stdout, stderr = proc.communicate()
    return stdout


def align_reference(path, index, sequences, nproc):
    cmd = ('bowtie2 --local -x {} -U {} -f --no-unal --no-sq --no-head '
           '--ma 5 -p {}').format(index, sequences, nproc)
    proc = subprocess.Popen(shlex.split(cmd),
                            stderr=subprocess.PIPE,
                            stdout=subprocess.PIPE,
                            cwd=path)

    stdout, stderr = proc.communicate()
    return stdout


def get_reader(output):
    fieldnames = [
        'seq_id',
        'flags',
        'reference',
        'ref_offset',
        'map_quality',
        'cigar',
        'mate_reference',
        'ref_mate_offset',
        'frag_size',
        'read_seq',
        'read_quality',
    ]
    return csv.DictReader(StringIO.StringIO(output), delimiter='\t',
                          fieldnames=fieldnames, restkey='optional')


def create_seqs(read_seq, ref_seq, cigar, ref_offset, **kwargs):
    ref_offset = max(0, ref_offset - 1)
    skip_ref = ref_seq[:ref_offset]
    ref_seq = ref_seq[ref_offset:]

    final_seq = []
    final_ref = []
    skips = []
    for (cnt, op) in re.findall('(\d+)([A-Z])', cigar):
        cnt = int(cnt)
        if op == 'M':
            final_seq.extend(read_seq[:cnt])
            final_ref.extend(ref_seq[:cnt])
            read_seq = read_seq[cnt:]
            ref_seq = ref_seq[cnt:]
        elif op == 'S':
            skips.append(read_seq[:cnt])
            read_seq = read_seq[cnt:]
        elif op == 'D':
            final_seq.extend(['-'] * cnt)
            final_ref.extend(ref_seq[:cnt])
            ref_seq = ref_seq[cnt:]
        elif op == 'I':
            final_ref.extend(['-'] * cnt)
            final_seq.extend(read_seq[:cnt])
            read_seq = read_seq[cnt:]
        else:
            raise Exception('Unknown opcode {}'.format(op))

    final_ref = skip_ref + ''.join(final_ref)
    final_seq = ''.join(final_seq)

    return (final_ref, final_seq, skips, skip_ref)


def add_imgt_gaps(imgt_germline, aligned_germline, sequence):
    for pos, size in gap_positions(imgt_germline):
        pos += gaps_before(gap_positions(aligned_germline), pos)
        aligned_germline = ''.join((
            aligned_germline[:pos],
            GAP_PLACEHOLDER * size,
            aligned_germline[pos:]
        ))
        sequence = ''.join((
            sequence[:pos],
            GAP_PLACEHOLDER * size,
            sequence[pos:]
        ))

    cdr3_start = 0
    count = 0
    for i, c in enumerate(aligned_germline):
        if c != '-':
            count += 1
            if count > CDR3_OFFSET:
                cdr3_start = i
                break
    else:
        cdr3_start = CDR3_OFFSET

    return (aligned_germline, sequence, cdr3_start)


def get_formatted_ties(genes):
    res = {}
    for ties, seq in genes.iteritems():
        res[format_ties(ties)] = seq
    return res


def process_sample(session, sample, indexes, temp, v_germlines, j_germlines,
                   nproc):
    indels = session.query(
        Sequence.ai,
        Sequence.seq_id,
        Sequence.sample_id,
        Sequence.sequence
    ).filter(
        Sequence.sample_id == sample.id,
        Sequence.probable_indel_or_misalign == 1
    )
    # Get the sequences that were not identifiable
    noresults = session.query(NoResult).filter(
        NoResult.sample_id == sample.id)

    if indels.count() == 0 and noresults.count() == 0:
        logger.info('Sample {} has no indels or noresults'.format(
            sample.id))
        return
    logger.info('Sample {} has {} indels and {} noresults'.format(
                sample.id, indels.count(), noresults.count()))

    mut_bucket = v_germlines.mut_bucket(sample.v_ties_mutations)
    len_bucket = v_germlines.length_bucket(sample.v_ties_len)
    bucket = '{}_{}'.format(str(mut_bucket).replace('.', ''),
                            len_bucket)
    sample_v_germlines = get_formatted_ties(v_germlines.all_ties(
            sample.v_ties_len, sample.v_ties_mutations))
    sample_j_germlines = get_formatted_ties(j_germlines.all_ties(
        sample.v_ties_len, sample.v_ties_mutations))
    if bucket not in indexes:
        indexes.add(bucket)
        v_path = os.path.join(temp, 'v_genes_{}'.format(bucket))
        j_path = os.path.join(temp, 'j_genes_{}'.format(bucket))
        logger.info('Creating index for V-ties at {} length, {} '
                    'mutation'.format(len_bucket, mut_bucket))
        build_index(sample_v_germlines, v_path)
        build_index(sample_j_germlines, j_path)

    seq_path = os.path.join(temp, 'll_{}.fasta'.format(sample.id))
    with open(seq_path, 'w+') as fh:
        fh.write(get_fasta({'tp=Sequence|ai={}|sample_id={}|seq_id={}'.format(
                r.ai, r.sample_id, r.seq_id): r.sequence for r in indels}))
        fh.write(get_fasta({'tp=NoResult|pk={}|sample_id={}|seq_id={}'.format(
            r.pk, r.sample_id, r.seq_id): r.sequence for r in noresults}))

    alignments = {}
    logger.info('Running bowtie2 for V-gene sequences')
    for line in get_reader(align_reference(temp, 'v_genes_{}'.format(bucket),
                                           seq_path, nproc)):
        line['ref_offset'] = int(line['ref_offset'])
        ref_gene = line['reference']
        ref, seq, rem_seqs, skip_ref = create_seqs(
            ref_seq=sample_v_germlines[ref_gene].replace('-', ''), **line)
        if len(rem_seqs) == 0:
            continue
        seq = ('N' * len(skip_ref)) + seq

        ref, seq, cdr3_start = add_imgt_gaps(sample_v_germlines[ref_gene], ref,
                                             seq)
        alignments[line['seq_id']] = {
            'v_germline': ref[:cdr3_start],
            'v_gene': line['reference'],
            'padding': len(skip_ref),
            'v_sequence': seq,
            'v_rem_seq': rem_seqs[-1],
            'cdr3_start': CDR3_OFFSET
        }

    seq_path = os.path.join(temp, 'll_j_{}.fasta'.format(sample.id))
    with open(seq_path, 'w+') as fh:
        fh.write(
            get_fasta({k: v['v_rem_seq'] for k, v in alignments.iteritems()}))

    tasks = []
    logger.info('Running bowtie2 for J-gene sequences')
    for line in get_reader(align_reference(temp, 'j_genes_{}'.format(bucket),
                                           seq_path, nproc)):
        line['ref_offset'] = int(line['ref_offset'])
        ref_gene = line['reference']
        ref, seq, rem_seqs, skip_ref = create_seqs(
            ref_seq=sample_j_germlines[ref_gene].replace('-', ''), **line)
        alignments[line['seq_id']]['j_gene'] = line['reference']

        full_seq = (alignments[line['seq_id']]['v_sequence'] +
                    alignments[line['seq_id']]['v_rem_seq'])

        cdr3_end = len(full_seq)
        for i in range(j_germlines.upstream_of_cdr3):
            if ref[-i] != '-':
                cdr3_end -= 1
        alignments[line['seq_id']]['cdr3_end'] = cdr3_end

        cdr3_length = cdr3_end - alignments[line['seq_id']]['cdr3_start']
        j_start = cdr3_start + cdr3_length

        full_germ = (alignments[line['seq_id']]['v_germline'] +
                     (GAP_PLACEHOLDER * cdr3_length))
        j_length = len(full_seq) - len(full_germ)
        if j_length <= 0:
            continue
        full_germ += ref[-j_length:]

        r_type, pk, sample_id, seq_id = [
            v.split('=', 1)[1] for v in line['seq_id'].split('|', 3)]
        insertions = gap_positions(full_germ)
        deletions = gap_positions(full_seq)

        alignment = VDJAlignment(
            VDJSequence(seq_id, full_seq.replace(GAP_PLACEHOLDER, '-'))
        )
        alignment.germline = full_germ.replace(GAP_PLACEHOLDER, '-')
        alignment.v_gene.add(GeneName(alignments[line['seq_id']]['v_gene']))
        alignment.j_gene.add(GeneName(alignments[line['seq_id']]['j_gene']))
        alignment.seq_offset = alignments[line['seq_id']]['padding']
        alignment.v_length = alignments[line['seq_id']]['cdr3_start']
        alignment.j_length = j_length
        alignment.v_mutation_fraction = 1 - (alignment.v_match /
                                             float(alignment.v_length))
        alignment.cdr3_start = alignments[line['seq_id']]['cdr3_start']
        alignment.cdr3_num_nts = cdr3_length
        alignment.post_cdr3_length = j_length
        alignment.insertions = insertions
        alignment.deletions = deletions
        alignment.locally_aligned = True

        tasks.append({
            'r_type': r_type,
            'pk': int(pk),
            'sample_id': int(sample_id),
            'alignment': alignment
        })
    return tasks


def add_sequences_from_sample(session, sample, sequences, props):
    for sequence in periodic_commit(session, sequences):
        alignment = sequence['alignment']
        if sequence['r_type'] == 'NoResult':
            try:
                props.validate(alignment)
            except AlignmentException as e:
                logger.info('NoResult could not be aligned {}'.format(
                    alignment.sequence.ids[0]))
                continue
            add_as_sequence(session, alignment, sample)
        elif sequence['r_type'] == 'Sequence':
            fields = {
                'partial': alignment.partial,

                'probable_indel_or_misalign': alignment.has_possible_indel,

                'v_gene': format_ties(alignment.v_gene),
                'j_gene': format_ties(alignment.j_gene),

                'num_gaps': alignment.num_gaps,
                'pad_length': alignment.pad_len,

                'v_match': alignment.v_match,
                'v_length': alignment.v_length,
                'j_match': alignment.j_match,
                'j_length': alignment.j_length,

                'removed_prefix': alignment.sequence.removed_prefix_sequence,
                'removed_prefix_qual':
                    alignment.sequence.removed_prefix_quality,
                'v_mutation_fraction': alignment.v_mutation_fraction,

                'pre_cdr3_length': alignment.pre_cdr3_length,
                'pre_cdr3_match': alignment.pre_cdr3_match,
                'post_cdr3_length': alignment.post_cdr3_length,
                'post_cdr3_match': alignment.post_cdr3_match,

                'in_frame': alignment.in_frame,
                'functional': alignment.functional,
                'stop': alignment.stop,
                'copy_number': len(alignment.sequence.ids),

                'cdr3_nt': alignment.cdr3,
                'cdr3_num_nts': len(alignment.cdr3),
                'cdr3_aa': lookups.aas_from_nts(alignment.cdr3),

                'sequence': str(alignment.sequence.sequence),
                'quality': alignment.sequence.quality,

                'locally_aligned': alignment.locally_aligned,
                '_insertions': serialize_gaps(alignment.insertions),
                '_deletions': serialize_gaps(alignment.deletions),

                'germline': alignment.germline
            }
            session.query(Sequence).filter(
                Sequence.ai == sequence['pk']
            ).update(fields, synchronize_session=False)


def remove_duplicates(session, sample):
    seqs = session.query(
        Sequence.ai, Sequence.seq_id, Sequence.v_gene, Sequence.j_gene,
        Sequence.cdr3_num_nts, Sequence.copy_number, Sequence.sequence
    ).filter(
        Sequence.locally_aligned.is_(True),
        Sequence.sample_id == sample.id
    ).order_by(Sequence.ai)

    for seq in seqs:
        potential_collapse = session.query(
            Sequence.ai, Sequence.sequence
        ).filter(
            Sequence.sample_id == sample.id,
            Sequence.v_gene == seq.v_gene,
            Sequence.j_gene == seq.j_gene,
            Sequence.cdr3_num_nts == seq.cdr3_num_nts,
        ).order_by(desc(Sequence.copy_number), Sequence.ai)

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
                    sample_id=sample.id
                ))
                session.query(Sequence).filter(Sequence.ai == seq.ai).delete()
                break
    session.commit()


def run_fix_sequences(session, args):
    v_germlines = VGermlines(args.v_germlines)
    j_germlines = JGermlines(args.j_germlines, args.upstream_of_cdr3)

    indexes = set()
    props = IdentificationProps(**args.__dict__)
    tasks = concurrent.TaskQueue()
    for sample in session.query(Sample):
        sequences = process_sample(session, sample, indexes, args.temp,
                                   v_germlines, j_germlines, args.nproc)
        add_sequences_from_sample(session, sample, sequences, props)
        remove_duplicates(session, sample)
