import argparse
import json

import numpy as np

from sqlalchemy import distinct, func

import sldb.common.config as config
from sldb.common.models import Clone, NoResult, SampleStats, Sequence
import sldb.util.lookups as lookups


_dist_fields = [
    'v_match',
    'j_match',
    'j_length',
    'v_gene',
    'j_gene',
    'copy_number',

    ('v_length', lambda seq: seq.v_length + seq.num_gaps),
    ('v_identity', lambda seq: np.ceil(100 * seq.v_match / seq.v_length)),
    ('cdr3_length', lambda seq: seq.junction_num_nts),
]

_seq_contexts = {
    'all': {
        'record_filter': lambda seq: True,
        'use_copy': True
    },
    'functional': {
        'record_filter': lambda seq: seq.functional,
        'use_copy': True
    },
    'nonfunctional': {
        'record_filter': lambda seq: not seq.functional,
        'use_copy': True
    },
    'unique': {
        'record_filter': lambda seq: seq.functional,
        'use_copy': False
    },
    'unique_multiple': {
        'record_filter': lambda seq: seq.functional and seq.copy_number > 1,
        'use_copy': False
    },
}

_clone_contexts = {
    'clones_all': {
        'record_filter': lambda functional: True,
    },
    'clones_functional': {
        'record_filter': lambda functional: functional
    },
    'clones_nonfunctional': {
        'record_filter': lambda functional: not functional
    },
}

class ContextStats(object):
    def __init__(self, record_filter):
        self._record_filter = record_filter
        self.distributions = {}

        for dist in _dist_fields:
            self.distributions[
                dist[0] if isinstance(dist, tuple) else dist] = {}

        self.sequence_cnt = 0
        self.in_frame_cnt = 0
        self.stop_cnt = 0
        self.functional_cnt = 0


class SeqContextStats(ContextStats):
    def __init__(self, record_filter, use_copy):
        super(SeqContextStats, self).__init__(record_filter)
        self._use_copy = use_copy

    def add_if_match(self, seq):
        if not self._record_filter(seq):
            return

        self.sequence_cnt += 1
        if seq.in_frame:
            self.in_frame_cnt += 1
        if seq.stop:
            self.stop_cnt += 1
        if seq.functional:
            self.functional_cnt += 1

        for dist in _dist_fields:
            if isinstance(dist, tuple):
                name = dist[0]
                value = dist[1](seq)
            else:
                name = dist
                value = getattr(seq, dist)

            if name not in self.distributions:
                self.distributions[name] = {}
            if value not in self.distributions[name]:
                self.distributions[name][value] = 0
            if self._use_copy:
                self.distributions[name][value] += seq.copy_number
            else:
                self.distributions[name][value] += 1


class CloneContextStats(ContextStats):
    def __init__(self, record_filter, seqs):
        super(CloneContextStats, self).__init__(record_filter)
        self._seqs = seqs

    def add_if_match(self, clone, in_frame, stop, functional):
        if not self._record_filter(functional):
            return

        self.sequence_cnt += 1
        if in_frame:
            self.in_frame_cnt += 1
        if stop:
            self.stop_cnt += 1
        if functional:
            self.functional_cnt += 1

        for dist in _dist_fields:
            name = dist[0] if isinstance(dist, tuple) else dist
            value = getattr(clone, name)
            if not isinstance(value, str):
                value = float(value)

            if name not in self.distributions:
                self.distributions[name] = {}
            if value not in self.distributions[name]:
                self.distributions[name][value] = 0
            self.distributions[name][value] += 1


def _get_cdr3_bounds(session, sample_id):
    cdr3_fld = Sequence.junction_num_nts
    cdr3s = []
    for seq in session.query(
                func.sum(Sequence.copy_number).label('copy_number'),
                cdr3_fld.label('cdr3_len')
            ).filter(
                Sequence.sample_id == sample_id,
                Sequence.copy_number > 1,
                Sequence.functional == 1
            ).group_by(cdr3_fld):
        cdr3s += [seq.cdr3_len] * int(seq.copy_number)
    if len(cdr3s) == 0:
        return None, None
    q25, q75 = np.percentile(cdr3s, [25, 75])
    iqr = q75 - q25
    return float(q25 - 1.5 * iqr), float(q75 + 1.5 * iqr)


def _add_stat(session, stats, sample_id, include_outliers,
              only_full_reads):
    for name, stat in stats.iteritems():
        ss = SampleStats(
            sample_id=sample_id,
            filter_type=name,
            outliers=include_outliers,
            full_reads=only_full_reads,
            sequence_cnt=stat.sequence_cnt,
            in_frame_cnt=stat.in_frame_cnt,
            stop_cnt=stat.stop_cnt,
            functional_cnt=stat.functional_cnt,
            no_result_cnt=session.query(NoResult).filter(
                NoResult.sample_id == sample_id).count()
        )
        for dname, dist in stat.distributions.iteritems():
            setattr(ss, '{}_dist'.format(dname),
                    json.dumps([(k, v) for k, v in dist.iteritems()]))
        session.add(ss)


def _calculate_seq_stats(session, sample_id, min_cdr3, max_cdr3,
                         include_outliers, only_full_reads):
    seq_statistics = {}
    for name, stat in _seq_contexts.iteritems():
        seq_statistics[name] = SeqContextStats(**stat)

    query = session.query(Sequence).filter(
        Sequence.sample_id == sample_id)
    if not include_outliers and min_cdr3 is not None:
        query = query.filter(Sequence.junction_num_nts >= min_cdr3,
                             Sequence.junction_num_nts <= max_cdr3)
    if only_full_reads:
        query = query.filter(Sequence.alignment == 'R1+R2')

    for seq in query:
        for stat in seq_statistics.values():
            stat.add_if_match(seq)

    _add_stat(session, seq_statistics, sample_id, include_outliers,
              only_full_reads)


def _calculate_clone_stats(session, sample_id, min_cdr3, max_cdr3,
                           include_outliers, only_full_reads):
    clone_statistics = {}
    for name, stat in _clone_contexts.iteritems():
        clone_statistics[name] = CloneContextStats(seqs=None, **stat)

    # TODO: This should be automatically generated from _dist_fields
    query = session.query(
        Sequence.clone_id,
        func.avg(Sequence.v_match).label('v_match'),
        func.avg(Sequence.j_match).label('j_match'),
        func.avg(Sequence.j_length).label('j_length'),
        Sequence.v_gene,
        Sequence.j_gene,
        func.avg(Sequence.copy_number).label('copy_number'),
        (Sequence.v_length + Sequence.num_gaps).label('v_length'),
        (
            func.ceil(100 * Sequence.v_match / Sequence.v_length)
        ).label('v_identity'),
        Sequence.junction_num_nts.label('cdr3_length')).filter(
            Sequence.sample_id == sample_id,
            ~Sequence.clone_id.is_(None)
        )

    if only_full_reads:
        query = query.filter(Sequence.alignment == 'R1+R2')
    query = query.group_by(Sequence.clone_id)
    for clone in query:
        clone_info = session.query(Clone.cdr3_nt).filter(
            Clone.id == clone.clone_id).first()
        in_frame = len(clone_info.cdr3_nt) % 3 == 0
        stop = '*' in lookups.aas_from_nts(clone_info.cdr3_nt, '')
        functional = in_frame and not stop
        for name, stat in clone_statistics.iteritems():
            stat.add_if_match(clone, in_frame, stop, functional)

    _add_stat(session, clone_statistics, sample_id, include_outliers,
              only_full_reads)


def _process_sample(session, sample_id, force):
    print 'Processing sample {}'.format(sample_id)
    existing_seq = session.query(Sequence).filter(
        Sequence.sample_id == sample_id)
    existing_nores = session.query(NoResult).filter(
        NoResult.sample_id == sample_id)
    if existing_seq.first() is None and existing_nores.first() is None:
        print '\tSKIPPING since there are no sequences in this DB version'
        return

    existing = session.query(SampleStats.sample_id).filter(
        SampleStats.sample_id == sample_id).first() is not None
    if force and existing:
        print '\tFORCING regeneration of stats'
        session.query(SampleStats).filter(
            SampleStats.sample_id == sample_id).delete()
        session.commit()

    if existing and not force:
        print ('\tSKIPPING stats since they already exists.'
               '  Use the --force flag to force regeneration.')
        return

    min_cdr3, max_cdr3 = _get_cdr3_bounds(session, sample_id)
    for include_outliers in [True, False]:
        print '\tOutliers={}'.format(include_outliers)
        for only_full_reads in [True, False]:
            print '\t\tOnly Full Reads={}'.format(only_full_reads)
            _calculate_seq_stats(session, sample_id, min_cdr3, max_cdr3,
                                 include_outliers, only_full_reads)
            _calculate_clone_stats(session, sample_id, min_cdr3, max_cdr3,
                                   include_outliers, only_full_reads)

            session.commit()


def run_sample_stats(session, args):
    if args.samples is None:
        samples = map(lambda s: s.id, session.query(Sample.id).all())
    else:
        samples = args.samples

    for sample_id in samples:
        _process_sample(session, sample_id, args.force)
