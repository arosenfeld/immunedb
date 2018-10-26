import itertools

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.spatial.distance as distance
import seaborn as sns

from immunedb.common.models import CloneStats, Sample
from immunedb.exporting.writer import ExportWriter
from immunedb.util.log import logger


def get_feature_str(sample, features):
    if 'sample' in features:
        return sample.name
    if 'subject' in features:
        return sample.subject.identifier
    s = ' '.join((sample.metadata_dict.get(f, 'NA') for f in features))
    return s


def get_sample_df(session, sample_ids, features, size_metric, dist_func):
    size_metric = {
        'copies': 'total_cnt',
        'instances': 'unique_cnt'
    }[size_metric]
    stats = session.query(
        CloneStats.sample_id,
        CloneStats.clone_id,
        getattr(CloneStats, size_metric).label('size')
    ).filter(
        CloneStats.sample_id.in_(sample_ids),
    )

    clone_sizes = {}
    for stat in stats:
        clone_sizes.setdefault(stat.clone_id, {})[stat.sample_id] = stat.size
    clone_sizes = pd.DataFrame(clone_sizes).fillna(0)

    sim = {}
    for s1, s2 in list(itertools.combinations(clone_sizes.index, 2)):
        fsim = 1 - dist_func(clone_sizes.loc[s1], clone_sizes.loc[s2])
        sim.setdefault(s1, {})[s2] = sim.setdefault(s2, {})[s1] = fsim

    return pd.DataFrame(sim)


def collapse_df_features(df, features, sample_instances, agg_func):
    def _apply_feature(arr):
        return [
            get_feature_str(sample_instances[sample_id], features)
            for sample_id in arr
        ]
    df.columns = _apply_feature(df.columns)
    df.index = _apply_feature(df.index)
    df = df.groupby(
        df.columns, axis=1
    ).agg(agg_func).groupby(
        df.index, axis=0
    ).agg(agg_func)
    return df.reindex(sorted(df.columns), axis=1).sort_index()


def write_clone_overlap(session, sample_ids=None, pool_on=('sample',),
                        size_metric='copies', sim_func='cosine',
                        agg_func='median', zipped=False, **kwargs):
    samples = session.query(Sample)
    if sample_ids:
        samples = samples.filter(Sample.id.in_(sample_ids))
    sample_instances = {s.id: s for s in samples}

    with ExportWriter(zipped=zipped) as writer:
        for subject in set([s.subject for s in sample_instances.values()]):
            logger.info('Calculating overlap for {}'.format(
                subject.identifier))
            sub_samples = [
                s.id for s in sample_instances.values() if s.subject == subject
            ]
            sdf = get_sample_df(session, sub_samples, pool_on, size_metric,
                                getattr(distance, sim_func))
            if sdf.empty:
                logger.warning(
                    'Subject {} had no clones for calculation'.format(
                        subject.identifier))
                continue

            sdf = collapse_df_features(sdf, pool_on, sample_instances,
                                       getattr(np, agg_func))
            name = '{}.overlap'.format(subject.identifier)

            with writer.get_handle(name + '.tsv') as fh:
                sdf.to_csv(fh, sep='\t')

            title_fmt = 'Subject {}\npooled by={}, similarity metric={}'
            if 'sample' not in pool_on:
                title_fmt += ', aggregation function={}'
            fig, ax = plt.subplots(figsize=(20, 20))
            ax = sns.heatmap(sdf, annot=True, linewidths=.25, vmin=0, vmax=1)
            ax.set_title(title_fmt.format(
                subject.identifier, ' & '.join(pool_on), sim_func, agg_func
            ))

            with writer.get_handle(name + '.pdf', 'wb+') as fh:
                plt.savefig(fh, bbox_inches='tight', format='pdf')

        return writer.get_zip_value()
