import csv

from sqlalchemy import distinct, func

from immunedb.common.models import Clone, CloneStats, Sequence, Subject
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
