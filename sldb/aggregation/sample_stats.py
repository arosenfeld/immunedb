import argparse
import json

import numpy as np

from sqlalchemy import distinct, func

import sldb.common.config as config
from sldb.common.models import *


_dist_fields = [
    Sequence.v_match,
    (Sequence.v_length + Sequence.num_gaps, 'v_length'),
    (func.ceil(100 * Sequence.v_match / Sequence.v_length), 'v_identity'),
    Sequence.j_match,
    Sequence.j_length,
    Sequence.v_gene,
    Sequence.j_gene,
    Sequence.copy_number,
    (Sequence.junction_num_nts, 'cdr3_length')
]

_seq_filters = [
    {
        'type': 'all',
        'filter_func': lambda q: q,
        'use_copy': True
    },
    {
        'type': 'functional',
        'filter_func': lambda q: q.filter(Sequence.functional == 1),
        'use_copy': True
    },
    {
        'type': 'nonfunctional',
        'filter_func': lambda q: q.filter(Sequence.functional == 0),
        'use_copy': True
    },
    {
        'type': 'unique',
        'filter_func': lambda q: q.filter(Sequence.functional == 1),
        'use_copy': False
    },
    {
        'type': 'unique_multiple',
        'filter_func': lambda q: q.filter(Sequence.functional == 1,
                                          Sequence.copy_number > 1),
        'use_copy': False
    },
]


_clone_filters = [
    {
        'type': 'clones_all',
        'filter_func': lambda q: q.filter(
            ~Sequence.clone_id.is_(None)),
        'use_copy': False,
    },
    {
        'type': 'clones_functional',
        'filter_func': lambda q: q.filter(
            ~Sequence.clone_id.is_(None),
            Sequence.functional == 1),
        'use_copy': False,
    },
    {
        'type': 'clones_nonfunctional',
        'filter_func': lambda q: q.filter(
            ~Sequence.clone_id.is_(None),
            Sequence.functional == 0),
        'use_copy': False,
    },
]


def _get_cdr3_bounds(session, filter_func, sample_id):
    cdr3_fld = Sequence.junction_num_nts
    cdr3s = []
    for seq in filter_func(session.query(
            func.sum(Sequence.copy_number).label('copy_number'),
            cdr3_fld.label('cdr3_len'))
            .filter(Sequence.sample_id == sample_id)
            .group_by(cdr3_fld)):
        cdr3s += [seq.cdr3_len] * int(seq.copy_number)
    if len(cdr3s) == 0:
        return None, None
    q25, q75 = np.percentile(cdr3s, [25, 75])
    iqr = q75 - q25
    return float(q25 - 1.5 * iqr), float(q75 + 1.5 * iqr)


def _get_clone_distribution(session, sample_id, key, filter_func, avg):
    result = {}
    if isinstance(key, tuple):
        key = key[0]
    q = filter_func(session.query(Sequence.seq_id,
                    func.avg(key).label('key') if avg else key.label('key')))\
        .filter(Sequence.sample_id == sample_id)\
        .group_by(Sequence.clone_id)

    for clone in q:
        key = float(clone.key) if avg else clone.key
        if key not in result:
            result[key] = 0
        result[key] += 1
    return result


def _get_distribution(session, sample_id, key, filter_func, use_copy,
                      cdr3_bounds, full_reads):
    result = {}
    if isinstance(key, tuple):
        key = key[0]

    q = filter_func(session.query(Sequence.seq_id,
                    key.label('key'))
                    .filter(Sequence.sample_id == sample_id))
    if full_reads:
        q = q.filter(Sequence.alignment == 'R1+R2')

    if use_copy:
        q = q.add_columns(func.sum(Sequence.copy_number).label('copy_number'))
    else:
        q = q.add_columns(func.count(Sequence.seq_id).label('unique'))

    if cdr3_bounds is not None:
        q = q.filter(Sequence.junction_num_nts >= cdr3_bounds[0],
                     Sequence.junction_num_nts <= cdr3_bounds[1])
    q = q.group_by(key)

    for row in q:
        result[row.key] = int(row.copy_number) if use_copy else \
            int(row.unique)
    return result


def _process_filter(session, sample_id, filter_type, filter_func,
                    use_copy, include_outliers, full_reads):
    if not include_outliers:
        min_cdr3, max_cdr3 = _get_cdr3_bounds(session, filter_func, sample_id)

    def base_query():
        q = filter_func(session.query(
            func.count(Sequence.seq_id).label('unique'))
            .filter(Sequence.sample_id == sample_id))
        if not include_outliers and min_cdr3 is not None:
            q = q.filter(Sequence.junction_num_nts >= min_cdr3,
                         Sequence.junction_num_nts <= max_cdr3)
        if full_reads:
            q = q.filter(Sequence.alignment == 'R1+R2')
        return q

    stat = SampleStats(sample_id=sample_id,
                       filter_type=filter_type,
                       outliers=include_outliers,
                       full_reads=full_reads)

    q = base_query().first()
    if q is None:
        stat.sequence_cnt = 0
    else:
        stat.sequence_cnt = q.unique

    q = base_query().filter(Sequence.in_frame == 1).first()
    if q is None:
        stat.in_frame_cnt = 0
    else:
        stat.in_frame_cnt = q.unique

    q = base_query().filter(Sequence.stop == 1).first()
    if q is None:
        stat.stop_cnt = 0
    else:
        stat.stop_cnt = q.unique

    q = base_query().filter(Sequence.functional == 1).first()
    if q is None:
        stat.functional_cnt = 0
    else:
        stat.functional_cnt = q.unique

    stat.no_result_cnt = (session.query(func.count(NoResult.seq_id))
                          .filter(NoResult.sample_id == sample_id).scalar())

    for dist in _dist_fields:
        dist_val = _get_distribution(
            session, sample_id, dist, filter_func, use_copy,
            ((min_cdr3, max_cdr3) if not include_outliers and min_cdr3
                is not None else None), full_reads)

        if isinstance(dist, tuple):
            name = dist[1]
        else:
            name = dist.name
        tuples = [[k, v] for k, v in sorted(dist_val.iteritems())]
        setattr(stat, '{}_dist'.format(name), json.dumps(tuples))
    session.add(stat)


def _process_clone_filter(session, sample_id, filter_type, filter_func,
                          include_outliers, full_reads):
    def base_query():
        q = filter_func(session.query(
            func.count(Sequence.seq_id).label('unique'))
            .filter(Sequence.sample_id == sample_id))\
            .group_by(Sequence.clone_id)
        if full_reads:
            q = q.filter(Sequence.alignment == 'R1+R2')
        return q

    stat = SampleStats(sample_id=sample_id,
                       filter_type=filter_type,
                       outliers=include_outliers,
                       full_reads=full_reads)
    stat.sequence_cnt = base_query().count()
    stat.in_frame_cnt = base_query().filter(Sequence.in_frame == 1).count()
    stat.stop_cnt = base_query().filter(Sequence.stop == 1).count()

    for dist in _dist_fields:
        avg = dist not in [Sequence.v_gene, Sequence.j_gene]
        dist_val = _get_clone_distribution(session, sample_id, dist,
                                           filter_func, avg)

        if isinstance(dist, tuple):
            name = dist[1]
        else:
            name = dist.name
        tuples = [[k, v] for k, v in sorted(dist_val.iteritems())]
        setattr(stat, '{}_dist'.format(name), json.dumps(tuples))
    session.add(stat)


def _process_sample(session, sample_id, force, clones_only):
    print 'Processing sample {}'.format(sample_id)
    existing_seq = session.query(Sequence).filter(
        Sequence.sample_id == sample_id)
    existing_nores = session.query(NoResult).filter(
        NoResult.sample_id == sample_id)
    if existing_seq.first() is None and existing_nores.first() is None:
        print '\tSKIPPING since there are no sequences in this DB version'
        return

    if force:
        print '\tFORCING regeneration of stats'
        if not clones_only:
            session.query(SampleStats).filter(
                SampleStats.sample_id == sample_id).delete()
        else:
            session.query(SampleStats).filter(
                SampleStats.sample_id == sample_id,
                SampleStats.filter_type.like('clones%')
            ).delete(synchronize_session='fetch')
        session.commit()

    existing = session.query(SampleStats).filter(
        SampleStats.sample_id == sample_id).first()
    if not force and existing is not None:
        print ('\tSKIPPING stats since they already exists.'
               '  Use the --force flag to force regeneration.')
        return

    for include_outliers in [True, False]:
        print '\tOutliers={}'.format(include_outliers)
        for full_reads in [True, False]:
            print '\t\tFull Reads={}'.format(full_reads)
            if not clones_only:
                for f in _seq_filters:
                    print ('\t\t\tGenerating sequence stats '
                           'for filter "{}"').format(
                        f['type'])
                    _process_filter(session, sample_id, f['type'],
                                    f['filter_func'], f['use_copy'],
                                    include_outliers, full_reads)
                    session.commit()

            for f in _clone_filters:
                print '\t\t\tGenerating sequence stats for filter "{}"'.format(
                    f['type'])
                _process_clone_filter(session, sample_id, f['type'],
                                      f['filter_func'], include_outliers,
                                      full_reads)
                session.commit()


def run_sample_stats(session, args):
    if args.samples is None:
        samples = map(lambda s: s.id, session.query(Sample.id).all())
    else:
        samples = args.samples

    for sample_id in samples:
        _process_sample(session, sample_id, args.force, args.clones_only)
