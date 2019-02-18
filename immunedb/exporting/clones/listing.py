from collections import Counter

import numpy as np
from sqlalchemy.orm import joinedload

from immunedb.common.models import CloneStats, Sample
from immunedb.exporting.tsv_writer import StreamingTSV
from immunedb.exporting.writer import ExportWriter
from immunedb.util.funcs import chunks
from immunedb.util.lookups import aas_from_nts
from immunedb.util.log import logger

DEFAULT_CLONE_FIELDS = [
    'clone_id', 'subject', 'v_gene', 'j_gene', 'functional', 'insertions',
    'deletions', 'cdr3_nt', 'cdr3_num_nts', 'cdr3_aa',
    'uniques', 'instances', 'copies', 'germline', 'parent_id',
    'avg_v_identity', 'top_copy_seq'
]


def get_clone_row(clone):
    row = {}
    for field in DEFAULT_CLONE_FIELDS:
        try:
            row[field] = getattr(clone, field)
        except AttributeError:
            pass
    row.update({
        'clone_id': clone.id,
        'subject': clone.subject.identifier,
        'functional': 'T' if clone.functional else 'F',
        'insertions': clone._insertions,
        'deletions': clone._deletions,
        'uniques': clone.overall_unique_cnt,
        'instances': clone.overall_instance_cnt,
        'copies': clone.overall_total_cnt,
    })
    return row


def get_immunedb_output(session, clones):
    writer = StreamingTSV(DEFAULT_CLONE_FIELDS)
    yield writer.writeheader()

    for clone, agg in clones.items():
        counts = agg['counts']
        row = get_clone_row(clone)
        row['copies'] = counts['copies']
        row['instances'] = counts['instances']
        row['top_copy_seq'] = agg['top_seq']
        row['avg_v_identity'] = round(agg['avg_v_identity'], 4)
        yield writer.writerow(row)


def get_vdjtools_output(session, clones):
    writer = StreamingTSV(['count', 'freq', 'cdr3nt', 'cdr3aa', 'v', 'd', 'j'])
    counts = Counter()
    total_copies = 0
    for clone, agg in clones.items():
        key = (clone.v_gene, clone.j_gene, clone.cdr3_nt)
        counts[key] += agg['counts']['copies']
        total_copies += counts[key]

    yield writer.writeheader()
    for key in sorted(counts, key=counts.get, reverse=True):
        count = counts[key]
        v, j, cdr3_nt = key
        yield writer.writerow({
            'count': count,
            'freq': count / total_copies,
            'cdr3nt': cdr3_nt,
            'cdr3aa': aas_from_nts(cdr3_nt),
            'v': v,
            'd': '.',
            'j': j,
        })


def _get_feature(stat, feature):
    if feature == 'subject':
        return stat.clone.subject.identifier
    elif feature == 'sample':
        return stat.sample.name
    else:
        return stat.sample.metadata_dict.get(feature, 'NA')


def get_pooled_samples(session, sample_ids, features):
    stats = session.query(
        CloneStats
    ).options(
        joinedload(CloneStats.sample),
        joinedload(CloneStats.clone).defer('tree'),
    )
    aggregated = {}
    for chunk_sample_ids in chunks(sample_ids, 20):
        for stat in stats.filter(CloneStats.sample_id.in_(chunk_sample_ids)):
            sample_feature = tuple(
                _get_feature(stat, f) for f in sorted(features)
            )
            key = (stat.sample.subject.identifier, sample_feature)
            agg = aggregated.setdefault(key, {}).setdefault(
                stat.clone,
                {'top_seq': None, 'top_copies': 0, 'counts': Counter(),
                    'avg_v_identity': []}
            )
            if stat.top_copy_seq_copies > agg['top_copies']:
                agg['top_copies'] = stat.top_copy_seq_copies
                agg['top_seq'] = stat.top_copy_seq_sequence
            agg['counts']['copies'] += stat.total_cnt
            agg['counts']['instances'] += stat.unique_cnt
            agg['avg_v_identity'].append(stat.avg_v_identity * stat.total_cnt)

    for clones in aggregated.values():
        for agg in clones.values():
            agg['avg_v_identity'] = np.sum(
                    agg['avg_v_identity']) / agg['counts']['copies']

    return aggregated


def get_filename(subject, feature_keys, feature_values):
    if feature_keys == ('subject',):
        return '{}.pooled.tsv'.format(subject)
    feature_value = '_AND_'.join(feature_values)
    return '{}.{}.pooled.tsv'.format(subject, feature_value)


def write_pooled_clones(session, out_format, sample_ids=None,
                        pool_on=('sample',), zipped=False, **kwargs):
    # Samples and subjects can't be combined with other features
    exclusives = set(pool_on).intersection(set(('sample', 'subject')))
    if len(pool_on) > 1 and exclusives:
        pool_on = (list(exclusives)[0],)
        logger.warning('You specified pooling on {feat} which '
                       'cannot be combined with other features.'
                       '  Using only {feat}.'.format(feat=pool_on[0]))

    logger.info('Writing clones pooled by {} in {} format'.format(
        ','.join(pool_on), out_format))

    sample_ids = sample_ids or [s.id for s in session.query(Sample)]
    aggregated = get_pooled_samples(session, sample_ids, pool_on)

    output_func = {
        'immunedb': get_immunedb_output,
        'vdjtools': get_vdjtools_output
    }[out_format]
    with ExportWriter(zipped=zipped) as fh:
        for (subject, feature_value), clones in aggregated.items():
            logger.info('Pooling subject {} for feature(s) {}'.format(
                subject,
                ','.join(feature_value)))
            fh.set_filename(get_filename(subject, pool_on, feature_value))
            fh.write(output_func(session, clones))
        return fh.get_zip_value()
