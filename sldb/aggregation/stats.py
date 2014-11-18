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
    (func.count(SequenceMapping.identity_seq_id), 'copy_number'),
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
        'filter_func': lambda q: q.filter(SequenceMapping.functional == 1)\
            .having(func.count(SequenceMapping.identity_seq_id) == 1),
        'use_copy': False
    },
    {
        'type': 'unique_multiple',
        'filter_func': lambda q: q.filter(
            SequenceMapping.functional == 1,
            func.count(SequenceMapping.identity_seq_id) > 1)
            .having(func.count(SequenceMapping.identity_seq_id) > 1),
        'use_copy': False
    },
]

_clone_filters = [
    {
        'type': 'clones_all',
        'filter_func': lambda q: q.filter(
            SequenceMapping.identity_seq.clone_id.isnot(None)),
        'use_copy': False,
    },
    {
        'type': 'clones_functional',
        'filter_func': lambda q: q.filter(
            SequenceMapping.identity_seq.clone_id.isnot(None),
            SequenceMapping.functional == 1),
        'use_copy': False,
    },
    {
        'type': 'clones_nonfunctional',
        'filter_func': lambda q: q.filter(
            SequenceMapping.identity_seq.clone_id.isnot(None),
            SequenceMapping.functional == 0),
        'use_copy': False,
    },
]


def _get_distribution(session, sample_id, key, filter_func, use_copy):
    if isinstance(key, tuple):
        key = key[0]

    result = {}
    q = filter_func(session.query(
            func.count(SequenceMapping).label('copy_number'),
            func.count(SequenceMapping).label('unique'))
        .join(Sequence)
        .filter(SequenceMapping.sample_id == sample_id))
    q = q.group_by(SequenceMapping.identity_seq_id)
        
    for row in q:
        result[row.key] = int(row.values)
    return result


def _process_filter(session, sample_id, filter_type, filter_func,
        use_copy):
    stat = SampleStats(sample_id=sample_id,
                       filter_type=filter_type)

    q = filter_func(session.query(
            func.count(SequenceMapping).label('copy_number'),
            func.count(SequenceMapping).label('unique'))
        .filter(SequenceMapping.sample_id == sample_id))
    q = q.group_by(SequenceMapping.identity_seq_id)

    if q is None:
        print 'Null stat for {}, {}'.format(sample_id, filter_type)
        return
    stat.sequence_cnt = q.first().copy_number if use_copy else q.first().unique
    '''
    stat.in_frame_cnt = filter_func(
        session.query(summation_func)\
            .filter(SequenceMapping.sample_id == sample_id,
                    SequenceMapping.in_frame)).scalar() or 0

    stat.stop_cnt = filter_func(
        session.query(summation_func)\
            .filter(SequenceMapping.sample_id == sample_id,
                    SequenceMapping.stop)).scalar() or 0
    for dist in _dist_fields:
        dist_val = _get_distribution(
            session, sample_id, dist, filter_func, use_copy)

        tuples = [[k, v] for k, v in sorted(dist_val.iteritems())]
        if isinstance(dist, tuple):
            dist = dist[1]
        setattr(stat, '{}_dist'.format(dist), json.dumps(tuples))
    '''
    session.add(stat)


def _process_sample(session, sample_id, force):
    print 'Processing sample {}'.format(sample_id)

    if force:
        print '\tFORCING regeneration of stats'
        session.query(SampleStats).filter(
            SampleStats.sample_id == sample_id).delete()
        session.query(CloneFrequency).filter(
            CloneFrequency.sample_id == sample_id).delete()
        session.commit()

    existing = session.query(SampleStats).filter(
        SampleStats.sample_id == sample_id).first()
    if not force and existing is not None:
        print ('\tSKIPPING stats since they already exists.'
               '  Use the --force flag to force regeneration.')
        return
    for f in _seq_filters + _clone_filters:
        print '\tGenerating sequence stats for filter "{}"'.format(f['type'])
        _process_filter(session, sample_id, f['type'], f['filter_func'],
                        f['use_copy'])

    for f in _clone_filters:
        for clone_info in f['filter_func'](session.query(
                Sequence.clone_id,
                func.count(Sequence.seq_id).label('unique'),
                func.sum(Sequence.copy_number_iden).label('total'))\
                .filter(Sequence.sample_id == sample_id,
                        Sequence.clone_id.isnot(None))\
                .group_by(Sequence.clone_id)):
            print '\tGenerating clone {} stats for filter "{}"'.format(
                clone_info.clone_id, f['type'])

            cf = CloneFrequency(
                sample_id=sample_id,
                clone_id=clone_info.clone_id,
                filter_type=f['type'],
                total_sequences=clone_info.total,
                unique_sequences=clone_info.unique)
            session.add(cf)
    session.commit()


def run_stats(session, args):
    if args.samples is None:
        samples = map(lambda s: s.id, session.query(Sample.id).all())
    else:
        samples = args.samples

    for sample_id in samples:
        _process_sample(session, sample_id, args.force)
