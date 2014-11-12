import argparse
import json

from sqlalchemy import func

import sldb.common.config as config
from sldb.common.models import *


_dist_fields = ['v_match', 'v_length', 'j_match', 'j_length', 'v_call',
                'j_call', 'v_gap_length', 'j_gap_length', 'copy_number_iden']


_seq_filters = [
    {
        'type': 'all',
        'filter_func': lambda q: q,
        'summation': True, 
    },
    {
        'type': 'functional',
        'filter_func': lambda q: q.filter(Sequence.functional),
        'summation': True, 
    },
    {
        'type': 'nonfunctional',
        'filter_func': lambda q: q.filter(not Sequence.functional),
        'summation': True, 
    },
    {
        'type': 'unique',
        'filter_func': lambda q: q.filter(not Sequence.functional,
                                          Sequence.copy_number_iden > 0),
        'summation': False, 
    },
    {
        'type': 'unique_multiple',
        'filter_func': lambda q: q.filter(not Sequence.functional,
                                          Sequence.copy_number_iden > 1),
        'summation': False, 
    },
]

class Stats(object):
    def __init__(self, session, sample_id):
        self.session = session
        self.sample_id = sample_id
        # Base statistics
        self.base_cnts = {k: 0 for k in filter(lambda e: e.endswith('_cnt'),
                                               dir(Sample))}
        # Functional specific statistics
        self.stats = [
            FilterStats(sample_id, 'all',
                        lambda e: True,
                        lambda e: e.copy_number_iden),
            FilterStats(sample_id, 'functional',
                        lambda e: e.functional,
                        lambda e: e.copy_number_iden),
            FilterStats(sample_id, 'nonfunctional',
                        lambda e: not e.functional,
                        lambda e: e.copy_number_iden),
            FilterStats(sample_id, 'unique',
                        lambda e: e.copy_number_iden > 0 and e.functional,
                        lambda e: 1),
            FilterStats(sample_id, 'unique_multiple',
                        lambda e: e.copy_number_iden > 1 and e.functional,
                        lambda e: 1)
        ]

        self.clone_stats = [
            FilterStats(sample_id, 'clones_all',
                        lambda e: True,
                        lambda e: 1),
            FilterStats(sample_id, 'clones_functional',
                        lambda e: e.functional,
                        lambda e: 1),
            FilterStats(sample_id, 'clones_nonfunctional',
                        lambda e: not e.functional,
                        lambda e: 1),
        ]

    def process_sequence(self, seq, size=None):
        self.base_cnts['valid_cnt'] += seq.copy_number_iden
        self.base_cnts['functional_cnt'] += \
            seq.copy_number_iden if seq.functional else 0

        for s in self.stats:
            if s.filter_fn(seq):
                s.process(seq)

        if seq.cluster is not None:
            for s in self.clone_stats:
                if s.filter_fn(seq):
                    s.process(seq)

                    cf = self.session.query(CloneFrequency).filter(
                        CloneFrequency.cluster == seq.cluster,
                        CloneFrequency.sample == seq.sample,
                        CloneFrequency.filter_type == s.filter_type).first()
                    if cf is not None:
                        cf.unique_sequences += 1
                        cf.total_sequences += seq.copy_number_iden
                    else:
                        cf = CloneFrequency(
                            sample=seq.sample,
                            cluster=seq.cluster,
                            filter_type=s.filter_type,
                            unique_sequences=1,
                            total_sequences=seq.copy_number_iden)

                    self.session.add(cf)

    def add_and_commit(self):
        self.session.query(Sample).filter_by(id=self.sample_id).update(
            self.base_cnts)
        for s in self.stats + self.clone_stats:
            self.session.add(s.get_populated_stat())
        self.session.commit()


class FilterStats(object):
    def __init__(self, sample_id, filter_type, filter_fn, count_fn):
        self._cnts = {}
        self._dists = {}
        self._stat = SampleStats(sample_id=sample_id, filter_type=filter_type)
        self.filter_type = filter_type
        self.filter_fn = filter_fn
        self.count_fn = count_fn

    def cnt_update(self, k, predicate, count):
        if k not in self._cnts:
            self._cnts[k] = 0
        if predicate:
            self._cnts[k] += count

    def dist_update(self, seq, stat):
        # Make sure the value of the statistic from the sequence has a count
        seq_val = getattr(seq, stat)
        if isinstance(seq_val, (int, long)):
            seq_val = int(seq_val)
        elif isinstance(seq_val, str):
            seq_val = seq_val.strip()
        elif seq_val is None:
            return

        self.dist_update_func(stat, seq_val, self.count_fn(seq))

    def dist_update_func(self, stat, value, count):
        if stat not in self._dists:
            self._dists[stat] = {}
        if value not in self._dists[stat]:
            self._dists[stat][value] = 0
        self._dists[stat][value] += count

    def get_populated_stat(self):
        for attrib in self._cnts:
            setattr(self._stat, attrib, self._cnts[attrib])
        for attrib in self._dists:
            tuples = [[k, v] for k, v in
                      sorted(self._dists[attrib].iteritems())]
            setattr(self._stat, '{}_dist'.format(attrib), json.dumps(tuples))
        return self._stat

    def process(self, e):
        self.cnt_update('in_frame_cnt', e.in_frame, self.count_fn(e))
        self.cnt_update('stop_cnt', e.stop, self.count_fn(e))
        self.cnt_update('sequence_cnt', True, self.count_fn(e))

        if e.junction_nt is not None:
            self.dist_update_func('cdr3_len', len(e.junction_nt),
                                  self.count_fn(e))

        for field in _dist_fields:
            self.dist_update(e, field)


def _get_distribution(session, sample_id, key, summation):
    if summation:
        agg = func.sum(Sequence.copy_number_iden).label('values')
    else:
        agg = func.count(Sequence.seq_id).label('values')

    result = {}
    for row in session.query(getattr(Sequence, key).label('key'), agg)\
            .filter(Sequence.sample_id == sample_id).group_by(key):
    result[row['key']] = row['value']
    return result


def _process_filter(session, sample_id, filter_type, filter_func, summation):
    stat = SampleStats(sample_id=sample_id,
                       filter_type=filter_type)

    for dist in _dist_fields:
        dist_val = _get_distribution(session, sample_id, dist, summation)

        tuples = [[k, v] for k, v in sorted(dist_val.iteritems())]
        setattr(stat, '{}_dist'.format(dist), json.dumps(tuples))
    return self._stat
    session.add(stat)


def _process_sample(session, sample_id):
    for f in _seq_filters:
        _process_filter(session, sample_id, f['type'], f['seq_filter_func'],
                        f['summation'])
    session.commit()


def run_stats():
    parser = config.get_base_arg_parser('Parse master-table into database.')
    parser.add_argument('--samples', nargs='+', type=int,
                        help='Limit statistics updates to certain samples')

    args = parser.parse_args()
    session = config.get_session(args)

    if args.samples is None:
        samples = map(lambda s: s.id, session.query(Sample.id).all())
    else:
        samples = args.samples

    for sample_id in samples:
        _process_sample(session, sample_id)
