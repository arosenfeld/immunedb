import csv

from sqlalchemy import func

from immunedb.common.models import Clone, CloneStats, Subject


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

    counts = session.query(
        CloneStats.clone_id, func.count(CloneStats.id).label('cnt')
    ).group_by(
        CloneStats.clone_id
    )
    counts = {c.clone_id: c.cnt - 1 for c in counts}
    sub_rows = {}
    for clone in session.query(Clone):
        sub_rows.setdefault(subjects[clone.subject_id], []).append({
            'count': clone.overall_unique_cnt,
            'cdr3nt': clone.cdr3_nt,
            'cdr3aa': clone.cdr3_aa,
            'v': clone.v_gene,
            'd': '.',
            'j': clone.j_gene,
            'VEnd': 0,
            'DStart': 0,
            'DEnd': 0,
            'JStart': 0,
            'incidence': counts[clone.id],
            'convergence': 1
        })

    for subject, rows in sub_rows.iteritems():
        total = sum([r['count'] for r in rows])
        rows = sorted(rows, key=lambda d: -d['count'])
        writers[subject].writeheader()
        for row in rows:
            row['freq'] = row['count'] / float(total)
            writers[subject].writerow(row)
