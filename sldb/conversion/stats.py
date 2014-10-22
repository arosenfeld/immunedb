import argparse
import json
from sldb.common.models import *


_dist_fields = ['v_match', 'v_length', 'j_match', 'j_length', 'v_call',
                'j_call', 'v_gap_length', 'j_gap_length', 'copy_number_close',
                'collapse_to_close', 'copy_number_iden', 'collapse_to_iden']


class Stats(object):
    def __init__(self, session, sample_id):
        self.session = session
        self.sample_id = sample_id
        ''' Base statistics '''
        self.base_cnts = {k: 0 for k in filter(lambda e: e.endswith('_cnt'),
                                               dir(Sample))}
        ''' Functional specific statistics '''
        self.stats = [
            FilterStats(sample_id, 'all', lambda e: True),
            FilterStats(sample_id, 'functional', lambda e:
                        e.functional),
            FilterStats(sample_id, 'nonfunctional', lambda e:
                        not e.functional),
            FilterStats(sample_id, 'unique', lambda e:
                        e.copy_number_iden > 0 and e.functional),
            FilterStats(sample_id, 'unique_multiple', lambda e:
                        e.copy_number_iden > 1 and e.functional),
        ]

        self.clone_stats = [
            FilterStats(sample_id, 'clones_all', lambda e:
                        e.clone_copy_number > 0),
            FilterStats(sample_id, 'clones_functional', lambda e:
                        e.clone_copy_number > 0 and e.functional),
            FilterStats(sample_id, 'clones_nonfunctional', lambda e:
                        e.clone_copy_number > 0 and not e.functional),
        ]

        self.clones_seen = set([])

    def process_sequence(self, seq, size=None, cn=None):
        self.base_cnts['valid_cnt'] += 1
        self.base_cnts['functional_cnt'] += 1 if seq.functional else 0

        for s in self.stats:
            if s.filter_fn(seq):
                s.process(seq)

        if seq.clone is not None and seq.clone not in self.clones_seen:
            for s in self.clone_stats:
                if s.filter_fn(seq):
                    s.process(seq)
                    cf = CloneFrequency(clone=seq.clone,
                                        sample_id=self.sample_id,
                                        size=size,
                                        copy_number=cn,
                                        filter_type=s.filter_type)
                    self.session.add(cf)
            self.clones_seen.add(seq.clone)

    def add_and_commit(self):
        self.session.query(Sample).filter_by(id=self.sample_id).update(
            self.base_cnts)
        for s in self.stats + self.clone_stats:
            self.session.add(s.get_populated_stat())
        self.session.commit()


class FilterStats(object):
    def __init__(self, sample_id, filter_type, filter_fn):
        self._cnts = {}
        self._dists = {}
        self._stat = SampleStats(sample_id=sample_id, filter_type=filter_type)
        self.filter_type = filter_type
        self.filter_fn = filter_fn

    def cnt_update(self, k, v):
        if k not in self._cnts:
            self._cnts[k] = 0
        if v:
            self._cnts[k] += 1

    def dist_update(self, seq, stat):
        # Make sure the value of the statistic from the sequence has a count
        seq_val = getattr(seq, stat)
        if isinstance(seq_val, (int, long)):
            seq_val = int(seq_val)
        elif isinstance(seq_val, str):
            seq_val = seq_val.strip()
        elif seq_val is None:
            return

        self.dist_set(stat, seq_val)

    def dist_set(self, stat, value):
        if stat not in self._dists:
            self._dists[stat] = {}
        if value not in self._dists[stat]:
            self._dists[stat][value] = 1
        else:
            self._dists[stat][value] += 1

    def get_populated_stat(self):
        for attrib in self._cnts:
            setattr(self._stat, attrib, self._cnts[attrib])
        for attrib in self._dists:
            tuples = [[k, v] for k, v in
                      sorted(self._dists[attrib].iteritems())]
            setattr(self._stat, '{}_dist'.format(attrib), json.dumps(tuples))
        return self._stat

    def process(self, e):
        self.cnt_update('in_frame_cnt', e.in_frame)
        self.cnt_update('stop_cnt', e.stop)
        self.cnt_update('mutation_inv_cnt', e.mutation_invariate)
        self.cnt_update('sequence_cnt', True)

        if e.junction_nt is not None:
            self.dist_set('cdr3_len', len(e.junction_nt))

        for field in _dist_fields:
            self.dist_update(e, field)
