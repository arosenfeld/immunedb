from collections import Counter, OrderedDict
import csv
import re
from sqlalchemy.orm import joinedload

from immunedb.identification.genes import GeneName
from immunedb.common.models import (Clone, CloneStats, Sample, Sequence,
                                    SequenceCollapse, Subject)
import immunedb.util.funcs as funcs
from immunedb.util.log import logger
from immunedb.util.lookups import aas_from_nts


def export_vdjtools(session, args):
    fieldnames = ['count', 'freq', 'cdr3nt', 'cdr3aa', 'v', 'd', 'j']
    if args.include_uniques:
        fieldnames.append('unique')

    clone_features = {c.id: (c.v_gene, c.j_gene, c.cdr3_nt)
                      for c in session.query(Clone.id, Clone.v_gene,
                                             Clone.j_gene, Clone.cdr3_nt)}
    for sample in session.query(Sample).order_by(Sample.id):
        logger.info('Exporting sample {}'.format(sample.name))
        sample_clones = {}
        stats = session.query(
            CloneStats.clone_id, CloneStats.total_cnt, CloneStats.unique_cnt
        ).filter(
            CloneStats.sample_id == sample.id
        )
        for stat in stats:
            key = clone_features[stat.clone_id]
            sample_clones.setdefault(key, Counter())['total'] += stat.total_cnt
            sample_clones[key]['unique'] += stat.unique_cnt

        writer = csv.DictWriter(
            open('{}.sample.txt'.format(sample.name), 'w+'),
            fieldnames=fieldnames,
            delimiter='\t',
            extrasaction='ignore'
        )
        total = float(sum([c['total'] for c in sample_clones.values()]))
        writer.writeheader()
        for key in sorted(sample_clones, key=sample_clones.get, reverse=True):
            counts = sample_clones[key]
            if counts['total'] < args.min_clone_size:
                continue
            v, j, cdr3_nt = key
            writer.writerow({
                'count': counts['total'],
                'freq': counts['total'] / total,
                'cdr3nt': cdr3_nt,
                'cdr3aa': aas_from_nts(cdr3_nt),
                'v': v,
                'd': '.',
                'j': j,
                'unique': counts['unique']
            })


def get_genbank_entries(seq, inference, gene_db):
    entry = OrderedDict()
    trim_start = re.search('[N-]*', seq.sequence).end()
    trim_end = re.search('[N-]*', seq.sequence[::-1]).end()
    trimmed_seq = seq.sequence[trim_start:len(seq.sequence) - trim_end]
    # Offsets of regions/segments
    v_gaps = trimmed_seq.count('-')
    trimmed_seq = trimmed_seq.replace('-', '')
    v_region = (1, len(trimmed_seq), 'V_region')
    v_segment = (1, 309 - v_gaps - trim_start, 'V_segment')
    cds = (
        v_segment[1] + 1,
        v_segment[1] + seq.cdr3_num_nts,
        'CDS'
    )
    j_segment = (cds[1] + 1, v_region[1], 'J_segment')

    # Key
    entry[('>Feature', seq.seq_id.split(' ')[0])] = []
    # V region
    entry[v_region] = []
    # V segment
    first_v = seq.v_gene.split('|')[0]
    entry[v_segment] = (
        ('gene', first_v),
        ('db_xref', '{}:{}'.format(gene_db, first_v)),
        ('inference', inference),
    )

    # J segment
    first_j = seq.j_gene.split('|')[0]
    entry[j_segment] = (
        ('gene', first_j),
        ('db_xref', '{}:{}'.format(gene_db, first_j)),
        ('inference', inference),
    )

    # CDS
    cds = ('<{}'.format(cds[0]), '>{}'.format(cds[1]), 'CDS')
    entry[cds] = (
        ('function', 'junction'),
        ('codon_start', 1),
        ('product', 'immunoglobulin heavy chain variable domain'),
        ('inference', inference)
    )

    return entry, trimmed_seq


def export_sample_genbank(session, sample_id, gene_db, inference, header):
    for sample in session.query(Sample):
        seqs = session.query(
            Sequence.seq_id,
            Sequence.v_gene,
            Sequence.j_gene,
            Sequence.cdr3_num_nts,
            Sequence.sequence,
            Sequence.copy_number
        ).filter(
            Sequence.sample_id == sample_id
        ).filter(Sequence.stop == 0)
        with open('{}.tbl'.format(sample.name), 'w+') as gb_fh:
            writer = csv.writer(gb_fh, delimiter='\t')
            with open('{}.fsa'.format(sample.name), 'w+') as fasta_fh:
                for seq in seqs:
                    gb_entry, fasta_seq = get_genbank_entries(
                        seq, inference, gene_db)
                    for entry, indented in gb_entry.iteritems():
                        writer.writerow(entry)
                        if indented:
                            for indent in indented:
                                writer.writerow(('', '', '') + indent)
                    seq_header = header + ' [note=AIRR_READ_COUNT:{}]'.format(
                        seq.copy_number)
                    fasta_fh.write('>{}\n{}\n'.format(
                        seq.seq_id.split(' ')[0] + ' {}'.format(seq_header),
                        fasta_seq))


def export_changeo(session, args):
    for subject in session.query(Subject):
        logger.info('Exporting subject {}'.format(subject.identifier))
        seqs = session.query(Sequence).filter(
            Sequence.subject_id == subject.id
        ).join(
            SequenceCollapse
        ).options(
            joinedload(Sequence.clone),
            joinedload(Sequence.collapse)
        )
        if args.clones_only:
            seqs = seqs.filter(~Sequence.clone_id.is_(None))
        if args.min_subject_copies is not None:
            seqs = seqs.filter(
                SequenceCollapse.copy_number_in_subject >=
                args.min_subject_copies
            )

        with open('{}.changeo.tsv'.format(subject.identifier), 'w+') as fh:
            writer = csv.DictWriter(
                fh, delimiter='\t', extrasaction='ignore', fieldnames=[
                    'SEQUENCE_ID',
                    'SEQUENCE_IMGT',
                    'FUNCTIONAL',
                    'IN_FRAME',
                    'STOP',
                    'V_CALL',
                    'J_CALL',
                    'JUNCTION_LENGTH',
                    'JUNCTION',
                    'V_SCORE',
                    'V_IDENTITY',
                    'J_SCORE',
                    'J_IDENTITY',
                    'DUPCOUNT',
                    'CLONE',
                    'CLONE_CDR3_NT'
                ]
            )
            writer.writeheader()
            for seq in funcs.yield_limit(seqs, Sequence.ai):
                v_prefix = GeneName(seq.v_gene).prefix
                j_prefix = GeneName(seq.j_gene).prefix
                writer.writerow({
                    'SEQUENCE_ID': seq.seq_id,
                    'SEQUENCE_IMGT': seq.sequence,
                    'FUNCTIONAL': 'T' if seq.functional else 'F',
                    'IN_FRAME': 'T' if seq.in_frame else 'F',
                    'STOP': 'T' if seq.stop else 'F',
                    'V_CALL': seq.v_gene.replace('|', ',' + v_prefix),
                    'J_CALL': seq.j_gene.replace('|', ',' + j_prefix),
                    'JUNCTION_LENGTH': seq.cdr3_num_nts,
                    'JUNCTION': seq.cdr3_nt,
                    'V_SCORE': round(
                        100 * seq.v_match / float(seq.v_length), 2),
                    'V_IDENTITY': seq.v_match,
                    'J_SCORE': round(
                        100 * seq.j_match / float(seq.j_length), 2),
                    'J_IDENTITY': seq.j_match,
                    'DUPCOUNT': seq.collapse.copy_number_in_subject,
                    'CLONE': seq.clone_id,
                    'CLONE_CDR3_NT': seq.clone.cdr3_nt if seq.clone else ''
                })


def export_genbank(session, args):
    args.inference = 'alignment:' + args.inference

    header = (
        '[organism={}] '
        '[moltype={}] '
        '[keywords=AIRR]').format(
        args.species, args.mol_type)
    samples = args.sample_ids or [s.id for s in session.query(Sample.id).all()]

    for sample_id in samples:
        logger.info('Exporting sample {}'.format(sample_id))
        export_sample_genbank(session, sample_id, args.gene_db, args.inference,
                              header)
