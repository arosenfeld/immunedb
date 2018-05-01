from collections import Counter
import csv

from immunedb.common.models import Clone, CloneStats, Sample
from immunedb.util.log import logger
from immunedb.util.lookups import aas_from_nts


def write_vdjtools(session, args):
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
