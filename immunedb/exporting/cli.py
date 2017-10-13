from collections import OrderedDict
import csv
import re

from sqlalchemy import distinct, func

from immunedb.common.models import Clone, CloneStats, Sample, Sequence, Subject
from immunedb.util.log import logger


def export_vdjtools(session, args):
    fieldnames = [
        'count', 'freq', 'cdr3nt', 'cdr3aa', 'v', 'd', 'j', 'VEnd', 'DStart',
        'DEnd', 'JStart', 'incidence', 'convergence'
    ]

    writers = {}
    subjects = {}
    for subject in session.query(Subject.id, Subject.identifier):
        writers[subject.identifier] = csv.DictWriter(
            open('pool.{}.summary.txt'.format(subject.identifier), 'w+'),
            fieldnames=fieldnames, delimiter='\t')
        subjects[subject.id] = subject.identifier

    subq_samples = session.query(
        Sequence.clone_id,
        func.count(distinct(Sequence.sample_id)).label('num_samples'),
        func.count(Sequence.seq_id).label('instances'),
    ).group_by(Sequence.clone_id).subquery('subq_samples')
    query = session.query(
        Clone.v_gene, Clone.j_gene, Clone.cdr3_nt, Clone.cdr3_aa,
        Clone.subject_id,
        subq_samples.c.num_samples,
        subq_samples.c.instances,
        func.sum(Clone.overall_total_cnt).label('num_reads'),
        func.sum(Clone.overall_unique_cnt).label('num_uniques'),
        func.count(distinct(Clone.id)).label('num_clones')
    ).filter(
        Clone.id == subq_samples.c.clone_id
    ).group_by(
        Clone.v_gene, Clone.j_gene, Clone.cdr3_nt, Clone.subject_id
    )
    sub_rows = {}
    for clone in query:
        size = int(clone.num_reads)
        if size < args.min_clone_size:
            continue
        sub_rows.setdefault(subjects[clone.subject_id], []).append({
            'count': size,
            'cdr3nt': clone.cdr3_nt,
            'cdr3aa': clone.cdr3_aa,
            'v': clone.v_gene,
            'd': '.',
            'j': clone.j_gene,
            'VEnd': 0,
            'DStart': 0,
            'DEnd': 0,
            'JStart': 0,
            'incidence': int(clone.num_samples),
            'convergence': int(clone.num_clones)
        })

    for subject, rows in sub_rows.iteritems():
        total = sum([r['count'] for r in rows])
        rows = sorted(rows, key=lambda d: -d['count'])
        writers[subject].writeheader()
        for row in rows:
            row['freq'] = row['count'] / float(total)
            writers[subject].writerow(row)


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
    first_j = seq.v_gene.split('|')[0]
    entry[j_segment] = (
        ('gene', first_v),
        ('db_xref', '{}:{}'.format(gene_db, first_v)),
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
    sample = session.query(Sample).filter(Sample.id == sample_id).one()
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


def export_genbank(session, args):
    inference = 'alignment:' + args.inference

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
