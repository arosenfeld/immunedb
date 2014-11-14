import argparse
import json

from sqlalchemy import distinct, func

import sldb.common.config as config
from sldb.common.models import *


_dist_fields = ['v_match', 'v_length', 'j_match', 'j_length', 'v_call',
                'j_call', 'v_gap_length', 'j_gap_length', 'copy_number_iden',
                ('junction_nt', func.length, 'cdr3_length')]


_seq_filters = [
    {
        'type': 'all',
        'filter_func': lambda q: q,
        'summation_func': func.sum(Sequence.copy_number_iden),
    },
    {
        'type': 'functional',
        'filter_func': lambda q: q.filter(Sequence.functional == 1),
        'summation_func': func.sum(Sequence.copy_number_iden),
    },
    {
        'type': 'nonfunctional',
        'filter_func': lambda q: q.filter(Sequence.functional == 0),
        'summation_func': func.sum(Sequence.copy_number_iden),
    },
    {
        'type': 'unique',
        'filter_func': lambda q: q.filter(Sequence.functional == 1,
                                          Sequence.copy_number_iden == 1),
        'summation_func': func.count(Sequence.seq_id),
    },
    {
        'type': 'unique_multiple',
        'filter_func': lambda q: q.filter(Sequence.functional == 1,
                                          Sequence.copy_number_iden > 1),
        'summation_func': func.count(Sequence.seq_id),
    },
]

_clone_filters = [
    {
        'type': 'clones_all',
        'filter_func': lambda q: q.filter(Sequence.clone_id.isnot(None)),
        'summation_func': func.count(Sequence.seq_id),
    },
    {
        'type': 'clones_functional',
        'filter_func': lambda q: q.filter(Sequence.clone_id.isnot(None),
                                          Sequence.functional == 1),
        'summation_func': func.count(Sequence.seq_id),
    },
    {
        'type': 'clones_nonfunctional',
        'filter_func': lambda q: q.filter(Sequence.clone_id.isnot(None),
                                          Sequence.functional == 0),
        'summation_func': func.count(Sequence.seq_id),
    },
]


def _get_distribution(session, sample_id, key, filter_func, summation_func):
    if isinstance(key, tuple):
        key = key[1](getattr(Sequence, key[0]))
    else:
        key = getattr(Sequence, key)

    result = {}
    for row in filter_func(session.query(key.label('key'),
                                         summation_func.label('values'))\
            .filter(Sequence.sample_id == sample_id).group_by(key)):
        result[row.key] = int(row.values)
    return result


def _process_filter(session, sample_id, filter_type, filter_func,
        summation_func):
    stat = SampleStats(sample_id=sample_id,
                       filter_type=filter_type)

    stat.sequence_cnt = filter_func(
        session.query(summation_func)\
           .filter(Sequence.sample_id == sample_id)).scalar() or 0

    stat.in_frame_cnt = filter_func(
        session.query(summation_func)\
            .filter(Sequence.sample_id == sample_id,
                    Sequence.in_frame)).scalar() or 0

    stat.stop_cnt = filter_func(
        session.query(summation_func)\
            .filter(Sequence.sample_id == sample_id,
                    Sequence.stop)).scalar() or 0

    for dist in _dist_fields:
        dist_val = _get_distribution(
            session, sample_id, dist, filter_func, summation_func)

        tuples = [[k, v] for k, v in sorted(dist_val.iteritems())]
        if isinstance(dist, tuple):
            dist = dist[2]
        setattr(stat, '{}_dist'.format(dist), json.dumps(tuples))
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
                        f['summation_func'])

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
