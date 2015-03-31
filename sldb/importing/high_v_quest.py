import csv
import distance
import re

import sldb.common.modification_log as mod_log
from sldb.common.models import (DuplicateSequence, NoResult, Sample, Study,
                                Subject, Sequence)
from sldb.identification.v_genes import VGene, VGermlines, get_common_seq
import sldb.util.funcs as funcs
import sldb.util.lookups as lookups
import sldb.identification.germlines as j_germlines


def _get_value(line, key):
    _HEADERS = {
        'seq_id': 'Sequence ID',
        'v_gene': 'V-GENE and allele',
        'j_gene': 'J-GENE and allele',
        'v_info': 'V-REGION identity nt',
        'j_info': 'J-REGION identity nt',
        'functional': 'Functionality',
        'in_frame': 'JUNCTION frame',
        'functional_comment': 'Functionality comment',
        'indel': 'V-REGION potential ins/del',
        'cdr3': 'JUNCTION',
        'sequence': 'Sequence',

        'v_region': 'V-REGION',
        'junction': 'JUNCTION',
        'gapped_sequence': 'V-D-J-REGION',
    }
    return line[_HEADERS.get(key)]


def _match(s1, s2):
    return len(s1), len(s1) - distance.hamming(s1, s2)


class ImportException(Exception):
    pass


def _setup_import(session, args):
    study, new = funcs.get_or_create(session, Study,
                                     name=args.study_name)
    if new:
        print 'Created new study "{}"'.format(study.name)

    subject, new = funcs.get_or_create(session, Subject,
                                       study=study,
                                       identifier=args.subject)
    if new:
        print 'Created new subject "{}"'.format(subject.identifier)

    sample, new = funcs.get_or_create(session, Sample,
                                      study=study,
                                      name=args.sample_name)
    if new:
        print 'Created new sample "{}"'.format(sample.name)
        sample.date = args.date
        sample.subject = subject
        sample.subset = args.subset
        sample.tissue = args.tissue
        sample.disease = args.disease
        sample.lab = args.lab
        sample.experimenter = args.experimenter
    elif sample.subject != subject:
            raise ImportException('Tried to use existing sample with '
                                  'same name and different subject')
    session.commit()

    return sample


def _split_info(line):
    return map(int, line.split(' ')[0].split('/'))


def _pad_replace(seq):
    return seq.upper().replace('.', '-')


def _read_summary(session, summary_reader, read_type, sample):
    print 'Reading Summary'
    reads = {}
    avg_mutations = 0
    avg_length = 0
    for i, line in enumerate(summary_reader):
        if _get_value(line, 'sequence') is None:
            continue
        if i > 0 and i % 1000 == 0:
            print '\tRead {}'.format(i)
        seq_id = _get_value(line, 'seq_id')

        v_gene = map(
            lambda g: g[0],
            re.findall('(\d+-\d+((-\d)+)?(\*\d+)?)',
                       _get_value(line, 'v_gene')))

        j_gene = map(
            lambda g: g[0],
            re.findall('(\d+)(\*\d+)?',
                       _get_value(line, 'j_gene')))

        try:
            j_match, j_length = _split_info(_get_value(line, 'j_info'))
            v_match, v_length = _split_info(_get_value(line, 'v_info'))
        except:
            session.add(NoResult(
                seq_id=seq_id,
                sample=sample))
            continue
        reads[seq_id] = Sequence(
            seq_id=seq_id,
            sample_id=sample.id,

            alignment=read_type,
            probable_indel_or_misalign=len(_get_value(line, 'indel')) > 0,

            v_gene=v_gene,
            j_gene=j_gene,

            v_match=v_match,
            v_length=v_length,
            j_match=j_match,
            j_length=j_length,

            in_frame=_get_value(line, 'in_frame') == 'in-frame',
            functional=_get_value(line, 'functional').startswith('productive'),
            stop='stop codons' in _get_value(line, 'functional_comment'),
            copy_number=1,

            gap_method='IGMT',
        )
        avg_length += reads[seq_id].v_length
        avg_mutations += 1 - reads[seq_id].v_match / float(
            reads[seq_id].v_length)

    return reads, avg_mutations / float(len(reads)), avg_length / float(len(
        reads))


def _read_gapped(session, reads, gapped_reader, germlines, sample, use_v_ties,
                 muts, lens):
    print 'Reading Gapped'
    for i, line in enumerate(gapped_reader):
        if i > 0 and i % 1000 == 0:
            print '\tCommitted {}'.format(i)
            session.commit()
        seq_id = _get_value(line, 'seq_id')
        sequence = _get_value(line, 'gapped_sequence')
        if seq_id not in reads or sequence is None:
            continue
        sequence = _pad_replace(sequence)

        v_region = _pad_replace(_get_value(line, 'v_region'))
        junction = _get_value(line, 'junction')
        read = reads[seq_id]

        pad = re.match('-*', v_region)
        if pad is None:
            read.pad_length = 0
        else:
            read.pad_length = pad.end()
        read.num_gaps = v_region[read.pad_length:].count('-')

        read.junction_num_nts = len(junction)
        if len(junction) == 0:
            continue
        read.junction_nt = junction.upper()
        read.junction_aa = lookups.aas_from_nts(read.junction_nt, '')

        read.sequence = 'N' * read.pad_length + sequence[read.pad_length:]

        existing = session.query(Sequence).filter(
            Sequence.sequence == read.sequence,
            Sequence.sample == sample).first()

        if existing is not None:
            existing.copy_number += 1
            session.add(DuplicateSequence(duplicate_seq_id=existing.seq_id,
                                          sample_id=sample.id,
                                          seq_id=seq_id))
            continue

        try:
            if use_v_ties:
                read.v_gene = germlines.get_ties(read.v_gene, lens, muts)
            vs = [germlines['IGHV{}'.format(v)].sequence
                  for v in read.v_gene][:VGene.CDR3_OFFSET]
            read.germline = get_common_seq(vs)
            read.germline += '-' * read.junction_num_nts
            read.germline += j_germlines.j['IGHJ{}'.format(
                read.j_gene[0])][-j_germlines.j_offset:]
            read.germline = read.germline[:len(read.sequence)]
        except:
            session.add(NoResult(
                seq_id=read.seq_id,
                sample=sample))
            continue
        read.v_gene = funcs.format_ties(read.v_gene, 'IGHV')
        read.j_gene = funcs.format_ties(read.j_gene, 'IGHJ')

        read.sequence_replaced = ''.join(
            [g if s == 'N' else s for s, g in zip(
                read.sequence, read.germline)])

        try:
            read.pre_cdr3_length, read.pre_cdr3_match = _match(
                sequence[:VGene.CDR3_OFFSET],
                read.germline[:VGene.CDR3_OFFSET]
            )
            read.pre_cdr3_length -= read.pad_length
            read.post_cdr3_length, read.post_cdr3_match = _match(
                sequence[
                    VGene.CDR3_OFFSET + read.junction_num_nts:],
                read.germline[
                    VGene.CDR3_OFFSET + read.junction_num_nts:]
            )
        except:
            session.add(NoResult(
                seq_id=read.seq_id,
                sample=sample))

        session.add(read)
    session.commit()


def run_high_v_quest_import(session, args):
    try:
        sample = _setup_import(session, args)
    except ImportException as ex:
        print '[ERROR] Unable setup import: \n\t{}'.format(
            ex.message)
        return
    with open(args.summary_file) as summary_fh:
        reads, muts, lens = _read_summary(
            session,
            csv.DictReader(summary_fh, delimiter='\t'),
            args.read_type,
            sample)

    germlines = VGermlines(args.v_germlines, include_prepadded=True)

    mod_log.make_mod('hvquest_import', session=session, commit=True,
                     info=vars(args))

    with open(args.gapped_nt_file) as gapped_fh:
        _read_gapped(
            session,
            reads,
            csv.DictReader(gapped_fh, delimiter='\t'),
            germlines,
            sample,
            args.v_ties,
            muts,
            lens)
