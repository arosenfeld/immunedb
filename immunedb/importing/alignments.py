import csv
import itertools
import os

import dnautils

from immunedb.identification import (add_noresults_for_vdj, add_sequences,
                                     AlignmentException)
from immunedb.identification.genes import GeneName
from immunedb.identification.vdj_sequence import VDJAlignment, VDJSequence
from immunedb.identification.metadata import (parse_metadata, REQUIRED_FIELDS,
                                              MetadataException)
from immunedb.common.models import Sample, SampleMetadata, Study, Subject
import immunedb.util.funcs as funcs
from immunedb.util.log import logger
from Bio import SeqIO

from immunedb.identification.identify import IdentificationProps


def raw_germlines(fn):
    with open(fn) as fh:
        return {r.description: str(r.seq) for r in SeqIO.parse(fh, 'fasta')}


def create_sample(session, metadata):
    study, new = funcs.get_or_create(
        session, Study, name=metadata['study_name'])

    if new:
        logger.info('Created new study "{}"'.format(study.name))
        session.commit()

    sample, new = funcs.get_or_create(
        session, Sample, name=metadata['sample_name'], study=study)
    if new:
        logger.info('Created new sample "{}"'.format(sample.name))
        for key, value in metadata.items():
            if key not in REQUIRED_FIELDS:
                session.add(SampleMetadata(
                    sample=sample,
                    key=key,
                    value=value
                ))

        subject, new = funcs.get_or_create(
            session, Subject, study=study,
            identifier=metadata['subject'])
        sample.subject = subject
        session.commit()
    else:
        logger.error(
            'Sample "{}" already exists'.format(metadata['sample_name']))
        return
    return sample


def collapse_duplicates(alignment, uniques):
    seq_len = len(alignment.sequence.sequence)
    for other in uniques.setdefault(seq_len, []):
        if dnautils.equal(other.sequence.sequence,
                          alignment.sequence.sequence):
            other.sequence.copy_number += alignment.sequence.copy_number
            break
    else:
        uniques[seq_len].append(alignment)


def parse_file(fh, sample, session, alignment_func, props, v_germlines,
               j_germlines):
    uniques = {}
    reader = csv.DictReader(fh, delimiter='\t')
    for i, line in enumerate(reader):
        try:
            alignment = alignment_func(line, v_germlines, j_germlines)
            if props.trim_to:
                alignment.trim_to(props.trim_to)
            collapse_duplicates(alignment, uniques)
        except AlignmentException as e:
            add_noresults_for_vdj(session, e.args[0], sample, str(e.args[1]))

    metrics = {'muts': [], 'lens': []}
    for unique in itertools.chain.from_iterable(uniques.values()):
        try:
            add_sequences(session, [unique], sample)
            metrics['lens'].append(unique.v_length)
            metrics['muts'].append(unique.v_mutation_fraction)
        except AlignmentException as e:
            add_noresults_for_vdj(session, unique.sequence, sample, str(e))

    if metrics:
        sample.v_ties_len = sum(metrics['lens']) / len(metrics['lens'])
        sample.v_ties_mutations = sum(metrics['muts']) / len(metrics['muts'])


def add_imgt_gaps(germline, sequence):
    seq_gaps = [i for i, c in enumerate(sequence) if c == '.']
    added = 0
    for i, c in enumerate(germline):
        if c == '-':
            gaps_before = len([j for j in seq_gaps if j < i])
            sequence.add_gap(i + gaps_before, '.')
            added += 1
    return sequence, added


def parse_airr(line, v_germlines, j_germlines):
    seq = VDJSequence(
        seq_id=line['sequence_id'].replace('reversed|', ''),
        sequence=line['sequence_alignment'],
        rev_comp=line['rev_comp'] == 'T',
    )
    if not all([line['v_call'], line['j_call'], line['junction_aa']]):
        raise AlignmentException(seq, 'Missing v_gene, j_gene, or junction_aa')

    seq.pad(int(line['v_germline_start']) - 1)
    try:
        v_germ_seq = funcs.consensus(
            [v_germlines[v] for v in line['v_call'].split(',')]
        )
    except KeyError:
        raise AlignmentException(
            seq,
            'V-gene {} not in germline database'.format(line['v_call'])
        )

    aligned_germ = ''.join([
        v_germ_seq.replace('-', '')[:int(line['v_germline_start']) - 1],
        line['germline_alignment']
    ])
    # Append the missing portion, if any, of the J to the germline
    j_germ_seq = funcs.consensus(
        [j_germlines[j] for j in line['j_call'].split(',')]
    )
    append_j = len(j_germ_seq) - int(line['j_germline_end'])
    if append_j:
        aligned_germ += j_germ_seq[-append_j:]
        seq.pad_right(append_j)

    aligned_seq, gaps_added = add_imgt_gaps(v_germ_seq, seq)
    aligned_germ = add_imgt_gaps(
        v_germ_seq, VDJSequence('', aligned_germ)
    )[0].sequence
    cdr3_start = int(line['cdr3_start']) - int(line['v_sequence_start'])
    # Push the start of the CDR3 based on number of IMGT gaps added.  Then add
    # 3 because IgBLAST's CDR3 excludes the preserved Cysteine
    cdr3_start += gaps_added - 3
    cdr3_start += aligned_seq.sequence[:cdr3_start].count('-')
    cdr3_start += int(line['v_germline_start']) - 1
    cdr3_end = cdr3_start + len(line['cdr3']) + 6
    cdr3_seq = aligned_seq.sequence[cdr3_start:cdr3_end]

    aligned_germ = ''.join([
        aligned_germ[:cdr3_start],
        '.' * (cdr3_end - cdr3_start),
        aligned_germ[cdr3_end:]
    ])
    aligned_seq = ''.join([
        aligned_seq.sequence[:cdr3_start],
        cdr3_seq,
        aligned_seq.sequence[cdr3_end:]
    ])
    alignment = funcs.create_proxy(VDJAlignment(
        VDJSequence(line['sequence_id'], aligned_seq.replace('.', '-'))
    ))
    alignment.germline = aligned_germ.replace('.', '-')
    alignment.v_gene = set([GeneName(c) for c in line['v_call'].split(',')])
    alignment.j_gene = set([GeneName(c) for c in line['j_call'].split(',')])
    alignment.locally_aligned = True
    alignment.seq_offset = int(line['v_germline_start']) - 1
    alignment.v_length = int(line['v_alignment_end'])
    alignment.j_length = (int(line['j_alignment_end']) -
                          int(line['j_alignment_start']))
    alignment.v_mutation_fraction = (100 - float(line['v_identity'])) / 100
    alignment.cdr3_start = cdr3_start
    alignment.cdr3_num_nts = len(cdr3_seq)
    # Skipping the germline_cdr3 field and instead populating its dependencies
    # via the proxy
    alignment.v_match = float(line['v_identity']) * alignment.v_length / 100
    alignment.j_match = float(line['j_identity']) * alignment.j_length / 100
    alignment.post_cdr3_length = len(alignment.sequence.sequence) - cdr3_end
    alignment.insertions = funcs.gap_positions(aligned_germ)
    alignment.deletions = funcs.gap_positions(aligned_seq)

    return alignment


def import_alignments(session, args):
    parse_funcs = {
        'airr': parse_airr
    }

    meta_fn = args.metadata if args.metadata else os.path.join(
        args.sample_dir, 'metadata.tsv')
    if not os.path.isfile(meta_fn):
        logger.error('Metadata file not found.')
        return
    with open(meta_fn, 'rU') as fh:
        try:
            metadata = parse_metadata(session, fh, args.warn_existing,
                                      args.warn_missing, args.sample_dir)
        except MetadataException as ex:
            logger.error(ex.message)
            return

    props = IdentificationProps(**args.__dict__)
    v_germlines = raw_germlines(args.v_germlines)
    j_germlines = raw_germlines(args.j_germlines)

    for sample_name in sorted(metadata.keys()):
        sample = create_sample(session, metadata[sample_name])
        if sample:
            path = os.path.join(
                args.sample_dir, metadata[sample_name]['file_name'])
            with open(path) as fh:
                parse_file(fh, sample, session, parse_funcs[args.format],
                           props, v_germlines, j_germlines)
