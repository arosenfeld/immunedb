import csv
import dnautils
import os
import re

from immunedb.common.models import (CDR3_OFFSET, Sample, SampleMetadata, Study,
                                    Subject)
from immunedb.identification import (add_noresults_for_vdj, add_sequences,
                                     AlignmentException, get_common_seq)

from immunedb.identification.metadata import (MetadataException,
                                              parse_metadata, REQUIRED_FIELDS)
from immunedb.identification.genes import JGermlines, VGermlines, GeneName
from immunedb.identification.identify import IdentificationProps
from immunedb.identification.vdj_sequence import VDJAlignment, VDJSequence
import immunedb.util.funcs as funcs
from immunedb.util.log import logger


def create_alignment(seq, line, v_germlines, j_germlines):
    alignment = VDJAlignment(seq)
    alignment.v_gene = set([GeneName(g) for g in line['V_CALL'].split(',')])
    alignment.j_gene = set([GeneName(g) for g in line['J_CALL'].split(',')])
    alignment.cdr3_num_nts = int(line['JUNCTION_LENGTH'])
    alignment.v_length = int(line['V_SEQ_LENGTH'])
    alignment.seq_offset = re.match(r'[-]*', alignment.sequence.sequence).end()

    # TODO: Calculate these
    alignment.v_length = CDR3_OFFSET - seq[:CDR3_OFFSET].count('-')
    alignment.j_length = j_germlines.upstream_of_cdr3

    germ_v = [v_germlines[g] for g in alignment.v_gene if g in v_germlines]
    germ_j = [j_germlines[g] for g in alignment.j_gene if g in j_germlines]
    if len(germ_v) == 0 or len(germ_j) == 0:
        raise AlignmentException('Missing germlines: V={} J={}'.format(
            ','.join([str(v) for v in alignment.v_gene]),
            ','.join([str(j) for j in alignment.j_gene])))

    germ_v = get_common_seq(germ_v)
    germ_j = get_common_seq(germ_j)

    alignment.germline_cdr3 = ''.join((
        germ_v,
        '-' * (len(alignment.sequence) - len(germ_v) - len(germ_j)),
        germ_j
    ))[CDR3_OFFSET:CDR3_OFFSET + alignment.cdr3_num_nts]

    alignment.germline = ''.join([
        germ_v[:CDR3_OFFSET],
        '-' * alignment.cdr3_num_nts,
        germ_j[-j_germlines.upstream_of_cdr3:]
    ])

    alignment.sequence.pad_right(
        len(alignment.germline) -
        len(alignment.sequence.sequence)
    )

    if len(alignment.germline) != len(alignment.sequence.sequence):
        raise AlignmentException('Sequence and germline differ in size')
    return alignment


def extract_adaptive_sequence(idx, line, v_germlines, j_germlines):
    def _format_gene(g):
        g = g.replace('V0', 'V').replace('J0', 'J').replace('-0', '-')
        g = g.replace('TCRB', 'TRB')
        if '*' not in g:
            return g + '*01'
        return g

    def _resolve_gene(line, gene, db):
        full = _format_gene(line[gene.lower() + 'GeneName'])
        family = _format_gene(line[gene.lower() + 'FamilyName'])
        try:
            return full, db[GeneName(full)]
        except KeyError:
            return family, db[GeneName(family)]
        raise AlignmentException('Invalid {} gene: {} / {}'.format(
            gene.upper(), full, family))

    if not line['aminoAcid']:
        raise AlignmentException('No amino-acids provided')

    v_gene, v_germ = _resolve_gene(line, 'V', v_germlines)
    j_gene, j_germ = _resolve_gene(line, 'J', j_germlines)

    v_end = int(line['vIndex'])
    v_region = line['nucleotide'][:v_end]
    v_region = list(('N' * (CDR3_OFFSET - v_end) + v_region))

    for i in range(len(v_germ)):
        if v_germ[i] == '-':
            v_region[i] = '-'
    v_region = ''.join(v_region)
    cdr3_region = line['nucleotide'][v_end:v_end + int(line['cdr3Length'])]
    j_region = line['nucleotide'][v_end + int(line['cdr3Length']):]
    j_region = j_region + ('N' * (j_germlines.upstream_of_cdr3 -
                           len(j_region)))
    imgt_sequence = v_region + cdr3_region + j_region
    try:
        counts = line['count (templates/reads)']
    except KeyError:
        counts = line['count (templates)']

    return {
        'SEQUENCE_ID': 'seq_{}'.format(idx),
        'SEQUENCE_IMGT': imgt_sequence,
        'V_CALL': v_gene,
        'J_CALL': j_gene,
        'JUNCTION_LENGTH': line['cdr3Length'],
        'V_SEQ_LENGTH': v_end,
        'DUPCOUNT': counts
    }


def read_file(session, fmt, handle, sample, v_germlines, j_germlines, props):
    reader = csv.DictReader(handle, delimiter='\t')
    uniques = {}

    for i, line in enumerate(reader):
        if fmt == 'adaptive':
            try:
                line = extract_adaptive_sequence(i, line, v_germlines,
                                                 j_germlines)
            except (AlignmentException, KeyError) as e:
                seq = VDJSequence('seq_{}'.format(i), '')
                add_noresults_for_vdj(session, seq, sample, str(e))
                continue
        seq = VDJSequence(line['SEQUENCE_ID'],
                          line['SEQUENCE_IMGT'].replace('.', '-'))
        if 'DUPCOUNT' in line:
            seq.copy_number = int(line['DUPCOUNT'])
        try:
            alignment = create_alignment(seq, line, v_germlines, j_germlines)
            for other in uniques.setdefault(
                    len(alignment.sequence.sequence), []):
                if dnautils.equal(other.sequence.sequence,
                                  alignment.sequence.sequence):
                    other.sequence.copy_number += (
                        alignment.sequence.copy_number)
                    break
            else:
                uniques[len(alignment.sequence.sequence)].append(alignment)
        except AlignmentException as e:
            add_noresults_for_vdj(session, seq, sample, str(e))

    uniques = [s for k in sorted(uniques.keys()) for s in uniques[k]]
    lens = []
    muts = []
    for unique in uniques:
        try:
            props.validate(unique)
            add_sequences(session, [unique], sample)
            lens.append(unique.v_length)
            muts.append(unique.v_mutation_fraction)
        except AlignmentException as e:
            add_noresults_for_vdj(session, seq, sample, str(e))

    if len(lens) > 0:
        sample.v_ties_len = sum(lens) / len(lens)
        sample.v_ties_mutations = sum(muts) / len(muts)

    session.commit()


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


def run_import(session, args):
    v_germlines = VGermlines(args.v_germlines)
    j_germlines = JGermlines(args.j_germlines, args.upstream_of_cdr3,
                             args.anchor_len)

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
    for sample_name in sorted(metadata.keys()):
        sample = create_sample(session, metadata[sample_name])
        if sample:
            path = os.path.join(
                args.sample_dir, metadata[sample_name]['file_name'])
            with open(path) as fh:
                read_file(session, args.format, fh, sample, v_germlines,
                          j_germlines, props)
