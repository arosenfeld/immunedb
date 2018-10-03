from collections import Counter
import io
import zipfile

from sqlalchemy.sql import func
from sqlalchemy.orm import joinedload, defer

from immunedb.common.models import CloneStats, Sample, Sequence
from immunedb.exporting.tsv_writer import StreamingTSV, write_tsv
from immunedb.exporting.clones import DEFAULT_CLONE_FIELDS, get_clone_row
from immunedb.util.funcs import chunks
from immunedb.util.log import logger


def _get_output(session, clones, fields, top_seqs={}):
    writer = StreamingTSV(fields)
    yield writer.writeheader()
    for clone, agg in clones.items():
        counts = agg['counts']
        row = get_clone_row(clone, fields,
                            counts['mutations'] / counts['copies'])
        row['copies'] = counts['copies']
        row['instances'] = counts['instances']
        row['top_copy_seq'] = session.query(Sequence.sequence).filter(
            Sequence.ai == agg['top_seq']
        ).one().sequence
        yield writer.writerow(row)


def get_pooled_samples(session, sample_ids, feature):
    def _get_feature(stat):
        if feature == 'subject':
            return stat.clone.subject.identifier
        elif feature == 'sample':
            return stat.sample.name
        else:
            return stat.sample.metadata_dict.get(feature, 'NA')

    stats = session.query(
        CloneStats
    ).options(
        joinedload(CloneStats.sample),
        joinedload(CloneStats.clone).defer('tree'),
    )
    aggregated = {}
    for chunk_sample_ids in chunks(sample_ids, 20):
        for stat in stats.filter(CloneStats.sample_id.in_(chunk_sample_ids)):
            sample_feature = _get_feature(stat)
            key = (stat.sample.subject.identifier, sample_feature)
            agg = aggregated.setdefault(key, {}).setdefault(
                stat.clone,
                {'top_seq': None, 'top_copies': 0, 'counts': Counter()}
            )
            if stat.top_copy_seq_copies > agg['top_copies']:
                agg['top_copies'] = stat.top_copy_seq_copies
                agg['top_seq'] = stat.top_copy_seq_ai
            agg['counts']['copies'] += stat.total_cnt
            agg['counts']['instances'] += stat.unique_cnt
            agg['counts']['mutations'] += stat.total_mutations()

    return aggregated


def get_filename(subject, feature_key, feature_value):
    if feature_key == 'subject':
        return '{}.pooled.txt'.format(subject)
    return '{}.{}.pooled.txt'.format(subject, feature_value)


def get_fields():
    fields = DEFAULT_CLONE_FIELDS[:]
    fields.append('top_copy_seq')
    fields.remove('uniques')
    return fields


def write_pooled_samples(session, args):
    if not args.sample_ids:
        args.sample_ids = [s.id for s in session.query(Sample)]

    logger.info('Pooling samples by {}'.format(args.feature))
    aggregated = get_pooled_samples(session, args.sample_ids, args.feature)
    for (subject, feature_value), clones in aggregated.items():
        write_tsv(get_filename(subject, args.feature, feature_value),
                  _get_output, session, clones, get_fields())


def get_pooled_samples_as_zip(session, sample_ids, feature):
    aggregated = get_pooled_samples(session, sample_ids, feature)
    out = io.BytesIO()
    with zipfile.ZipFile(out, 'w') as fh:
        for (subject, feature_value), clones in aggregated.items():
            output = ''.join(_get_output(session, clones, get_fields()))
            fh.writestr(get_filename(subject, feature, feature_value), output)
    yield out.getvalue()
