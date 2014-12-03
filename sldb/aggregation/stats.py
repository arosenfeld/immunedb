import argparse
import json

from sqlalchemy import distinct, func

import sldb.common.config as config
from sldb.common.models import *


_dist_fields = [
    SequenceMapping.v_match,
    SequenceMapping.v_length,
    SequenceMapping.j_match,
    SequenceMapping.j_length,
    Sequence.v_call,
    Sequence.j_call,
    SequenceMapping.copy_number,
    (func.length(Sequence.junction_nt), 'cdr3_length')
]


_seq_filters = [
    {
        'type': 'all',
        'filter_func': lambda q: q,
        'use_copy': True
    },
    {
        'type': 'functional',
        'filter_func': lambda q: q.filter(SequenceMapping.functional == 1),
        'use_copy': True
    },
    {
        'type': 'nonfunctional',
        'filter_func': lambda q: q.filter(SequenceMapping.functional == 0),
        'use_copy': True
    },
    {
        'type': 'unique',
        'filter_func': lambda q: q.filter(SequenceMapping.functional == 1),
        'use_copy': False
    },
    {
        'type': 'unique_multiple',
        'filter_func': lambda q: q.filter(SequenceMapping.functional == 1,
                                          SequenceMapping.copy_number > 1),
        'use_copy': False
    },
]


_clone_filters = [
    {
        'type': 'clones_all',
        'filter_func': lambda q: q.filter(
            ~SequenceMapping.clone_id.is_(None)),
        'use_copy': False,
    },
    {
        'type': 'clones_functional',
        'filter_func': lambda q: q.filter(
            ~SequenceMapping.clone_id.is_(None),
            SequenceMapping.functional == 1),
        'use_copy': False,
    },
    {
        'type': 'clones_nonfunctional',
        'filter_func': lambda q: q.filter(
            ~SequenceMapping.clone_id.is_(None),
            SequenceMapping.functional == 0),
        'use_copy': False,
    },
]


def _get_clone_distribution(session, sample_id, key, filter_func, avg):
    result = {}
    if isinstance(key, tuple):
        key = key[0]
    q = filter_func(session.query(SequenceMapping.identity_seq_id,
            func.avg(key).label('key') if avg else key.label('key')))\
        .filter(SequenceMapping.sample_id == sample_id)\
        .join(Sequence)\
        .group_by(SequenceMapping.clone_id)
        
    for clone in q:
        key = float(clone.key) if avg else clone.key
        if key not in result:
            result[key] =0
        result[key] += 1
    return result


def _get_distribution(session, sample_id, key, filter_func, use_copy):
    result = {}
    if isinstance(key, tuple):
        key = key[0]
    q = filter_func(session.query(SequenceMapping.identity_seq_id,
            key.label('key'),
            func.count(SequenceMapping.identity_seq_id).label('unique'),
            func.sum(SequenceMapping.copy_number).label('copy_number'))
        .filter(SequenceMapping.sample_id == sample_id))\
        .join(Sequence)\
        .group_by(key)
        
    for row in q:
        result[row.key] = int(row.copy_number) if use_copy else \
            int(row.unique)
    return result


def _process_filter(session, sample_id, filter_type, filter_func,
        use_copy):
    def base_query():
        return filter_func(session.query(
            func.count(SequenceMapping.seq_id).label('unique'),
            func.sum(SequenceMapping.copy_number).label('copy_number'))\
            .filter(SequenceMapping.sample_id == sample_id))

    stat = SampleStats(sample_id=sample_id,
                       filter_type=filter_type)

    q = base_query().first()
    if q is None:
        stat.sequence_cnt = 0
    else:
        stat.sequence_cnt = q.copy_number if use_copy else q.unique

    q = base_query().filter(SequenceMapping.in_frame == 1).first()
    if q is None:
        stat.in_frame_cnt = 0
    else:
        stat.in_frame_cnt = q.copy_number if use_copy else q.unique

    q = base_query().filter(SequenceMapping.stop == 1).first()
    if q is None:
        stat.stop_cnt = 0
    else:
        stat.stop_cnt = q.copy_number if use_copy else q.unique

    for dist in _dist_fields:
        dist_val = _get_distribution(
            session, sample_id, dist, filter_func, use_copy)

        if isinstance(dist, tuple):
            name = dist[1]
        else:
            name = dist.name
        tuples = [[k, v] for k, v in sorted(dist_val.iteritems())]
        setattr(stat, '{}_dist'.format(name), json.dumps(tuples))
    session.add(stat)

def _process_clone_filter(session, sample_id, filter_type, filter_func):
    def base_query():
        return filter_func(session.query(
            func.count(SequenceMapping.seq_id).label('unique'))
            .filter(SequenceMapping.sample_id == sample_id))\
            .group_by(SequenceMapping.clone_id)

    stat = SampleStats(sample_id=sample_id,
                       filter_type=filter_type)
    stat.sequence_cnt = len(base_query().all())
    stat.in_frame_cnt = len(base_query()\
        .filter(SequenceMapping.in_frame == 1).all())
    stat.stop_cnt = len(base_query()\
        .filter(SequenceMapping.stop == 1).all())

    for dist in _dist_fields:
        avg = dist not in [Sequence.v_call, Sequence.j_call]
        dist_val = _get_clone_distribution(session, sample_id, dist,
                                           filter_func, avg)

        if isinstance(dist, tuple):
            name = dist[1]
        else:
            name = dist.name
        tuples = [[k, v] for k, v in sorted(dist_val.iteritems())]
        setattr(stat, '{}_dist'.format(name), json.dumps(tuples))
    session.add(stat)

def _process_sample(session, sample_id, force):
    print 'Processing sample {}'.format(sample_id)

    if force:
        print '\tFORCING regeneration of stats'
        session.query(SampleStats).filter(
            SampleStats.sample_id == sample_id).delete()
        session.commit()

    existing = session.query(SampleStats).filter(
        SampleStats.sample_id == sample_id).first()
    if not force and existing is not None:
        print ('\tSKIPPING stats since they already exists.'
               '  Use the --force flag to force regeneration.')
        return


    sample = session.query(Sample).filter(Sample.id == sample_id).first()
    sample.valid_cnt = session.query(func.count(SequenceMapping.seq_id))\
        .filter(SequenceMapping.sample_id == sample_id).scalar()
    sample.functional_cnt = session.query(func.count(SequenceMapping.seq_id))\
        .filter(SequenceMapping.sample_id == sample_id,
                SequenceMapping.functional == 1).scalar()
    sample.no_result_cnt = session.query(func.count(NoResult.seq_id))\
        .filter(NoResult.sample_id == sample_id).scalar()
    session.commit()

    for f in _seq_filters:
        print '\tGenerating sequence stats for filter "{}"'.format(f['type'])
        _process_filter(session, sample_id, f['type'], f['filter_func'],
                        f['use_copy'])
        session.commit()

    for f in _clone_filters:
        print '\tGenerating sequence stats for filter "{}"'.format(f['type'])
        _process_clone_filter(session, sample_id, f['type'], f['filter_func'])
        session.commit()


def run_stats(session, args):
    if args.samples is None:
        samples = map(lambda s: s.id, session.query(Sample.id).all())
    else:
        samples = args.samples

    for sample_id in samples:
        _process_sample(session, sample_id, args.force)
