import json

import numpy as np

from sqlalchemy import and_, func

import immunedb.common.config as config
import immunedb.common.modification_log as mod_log
from immunedb.common.models import (Clone, NoResult, Sample, SampleStats,
                                    SelectionPressure, Sequence)
import immunedb.util.concurrent as concurrent
import immunedb.util.lookups as lookups
from immunedb.util.log import logger


_dist_fields = [
    'v_match',
    'j_match',
    'j_length',
    'v_gene',
    'j_gene',
    'copy_number',
    'v_length',
    'v_identity',
    'cdr3_length',
    'sp_fwr',
    'sp_cdr'
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

    def _update(self, seq_record, in_frame, stop, functional, add):
        self.sequence_cnt += add
        if in_frame:
            self.in_frame_cnt += add
        if stop:
            self.stop_cnt += add
        if functional:
            self.functional_cnt += add

        for field in _dist_fields:
            try:
                value = getattr(seq_record, field)
                if value is None:
                    continue
            except AttributeError:
                continue
            if not isinstance(value, str):
                value = float(value)

            if field not in self.distributions:
                self.distributions[field] = {}
            if value not in self.distributions[field]:
                self.distributions[field][value] = 0

            self.distributions[field][value] += add


class SeqContextStats(ContextStats):
    def __init__(self, session, record_filter, use_copy):
        super(SeqContextStats, self).__init__(record_filter)
        self._session = session
        self._use_copy = use_copy
        self.quality = []

    def _update_quality(self, seq_record):
        if seq_record.quality is not None:
            diff = len(seq_record.quality) - len(self.quality)
            if diff > 0:
                self.quality.extend([[] for _ in range(0, diff)])

            for i, b in enumerate(seq_record.quality):
                if b != ' ':
                    self.quality[i].append(ord(b) - 33)

    def add_if_match(self, seq_record):
        if not self._record_filter(seq_record):
            return

        add = int(seq_record.copy_number) if self._use_copy else 1

        self._update(seq_record, seq_record.in_frame, seq_record.stop,
                     seq_record.functional, add)
        self._update_quality(seq_record)


class CloneContextStats(ContextStats):
    def __init__(self, record_filter, seqs):
        super(CloneContextStats, self).__init__(record_filter)
        self._seqs = seqs

    def add_if_match(self, clone, in_frame, stop, functional):
        if not self._record_filter(functional):
            return

        self._update(clone, in_frame, stop, functional, 1)


class SampleStatsWorker(concurrent.Worker):
    def __init__(self, session):
        self._session = session

    def do_task(self, args):
        if args['func'] == 'seq':
            func = self._calculate_seq_stats
        elif args['func'] == 'clone':
            func = self._calculate_clone_stats
        self.info(('Processing {} for sample {}, include_outliers {}, '
                   'only_full_reads {}').format(
                   'sequences' if args['func'] == 'seq' else 'clones',
                   args['sample_id'], args['include_outliers'],
                   args['only_full_reads']))
        func(args['sample_id'], args['min_cdr3'], args['max_cdr3'],
             args['include_outliers'], args['only_full_reads'])
        self._session.commit()

    def _add_stat(self, stats, sample_id, include_outliers,
                  only_full_reads):
        for name, stat in stats.items():
            ss = SampleStats(
                sample_id=sample_id,
                filter_type=name,
                outliers=include_outliers,
                full_reads=only_full_reads,
                sequence_cnt=stat.sequence_cnt,
                in_frame_cnt=stat.in_frame_cnt,
                stop_cnt=stat.stop_cnt,
                functional_cnt=stat.functional_cnt,
                no_result_cnt=self._session.query(NoResult).filter(
                    NoResult.sample_id == sample_id).count(),
            )

            if hasattr(stat, 'quality'):
                quality_dist = []
                for pos, quals in enumerate(stat.quality):
                    if len(quals) > 0:
                        quality_dist.append((pos, round(np.mean(quals), 2)))
                ss.quality_dist = json.dumps(quality_dist)
            else:
                ss.quality_dist = json.dumps([])

            for dname, dist in stat.distributions.items():
                setattr(ss, '{}_dist'.format(dname),
                        json.dumps([(k, v) for k, v in dist.items()]))
            self._session.add(ss)

    def _calculate_seq_stats(self, sample_id, min_cdr3, max_cdr3,
                             include_outliers, only_full_reads):
        seq_statistics = {}
        for name, stat in _seq_contexts.items():
            seq_statistics[name] = SeqContextStats(self._session, **stat)

        # TODO: This should be automatically generated from _dist_fields
        query = self._session.query(
            Sequence.quality,
            Sequence.sample_id,
            Sequence.v_match,
            Sequence.j_match,
            Sequence.j_length,
            Sequence.v_gene,
            Sequence.j_gene,
            Sequence.in_frame,
            Sequence.stop,
            Sequence.functional,
            Sequence.copy_number,
            (Sequence.v_length + Sequence.num_gaps).label('v_length'),
            (
                func.ceil(100 * Sequence.v_match / Sequence.v_length)
            ).label('v_identity'),
            Sequence.cdr3_num_nts.label('cdr3_length'),
        ).filter(
            Sequence.sample_id == sample_id
        )

        if not include_outliers and min_cdr3 is not None:
            query = query.filter(Sequence.cdr3_num_nts >= min_cdr3,
                                 Sequence.cdr3_num_nts <= max_cdr3)
        if only_full_reads:
            query = query.filter(Sequence.partial == 0)

        for seq in query:
            for name, stat in seq_statistics.items():
                stat.add_if_match(seq)

        self._add_stat(seq_statistics, sample_id, include_outliers,
                       only_full_reads)

    def _calculate_clone_stats(self, sample_id, min_cdr3, max_cdr3,
                               include_outliers, only_full_reads):
        clone_statistics = {}
        for name, stat in _clone_contexts.items():
            clone_statistics[name] = CloneContextStats(seqs=None, **stat)

        # TODO: This should be automatically generated from _dist_fields
        query = self._session.query(
            Sequence.clone_id,
            func.round(func.avg(Sequence.v_match)).label('v_match'),
            func.round(func.avg(Sequence.j_match)).label('j_match'),
            func.round(func.avg(Sequence.j_length)).label('j_length'),
            Sequence.v_gene,
            Sequence.j_gene,
            func.count(Sequence.seq_id).label('copy_number'),
            func.round(
                func.avg(Sequence.v_length + Sequence.num_gaps)
            ).label('v_length'),
            func.round(
                func.avg(100 * Sequence.v_match / Sequence.v_length)
            ).label('v_identity'),
            Sequence.cdr3_num_nts.label('cdr3_length'),
            SelectionPressure.sigma_fwr.label('sp_fwr'),
            SelectionPressure.sigma_cdr.label('sp_cdr'),
        ).join(
            SelectionPressure,
            and_(
                SelectionPressure.clone_id == Sequence.clone_id,
                SelectionPressure.sample_id == Sequence.sample_id
            ),
            isouter=True
        ).filter(
            Sequence.sample_id == sample_id,
            ~Sequence.clone_id.is_(None)
        )

        if only_full_reads:
            query = query.filter(Sequence.partial == 0)
        query = query.group_by(Sequence.clone_id)

        for clone in query:
            clone_info = self._session.query(Clone.cdr3_nt).filter(
                Clone.id == clone.clone_id).first()
            in_frame = len(clone_info.cdr3_nt) % 3 == 0
            stop = '*' in lookups.aas_from_nts(clone_info.cdr3_nt)
            functional = in_frame and not stop
            for name, stat in clone_statistics.items():
                stat.add_if_match(clone, in_frame, stop, functional)

        self._add_stat(clone_statistics, sample_id, include_outliers,
                       only_full_reads)


def _get_cdr3_bounds(session, sample_id):
    cdr3_fld = Sequence.cdr3_num_nts
    cdr3s = []

    query = session.query(
        func.sum(Sequence.copy_number).label('copy_number'),
        cdr3_fld.label('cdr3_len')
    ).filter(
        Sequence.sample_id == sample_id,
        Sequence.copy_number > 1,
        Sequence.functional == 1
    ).group_by(cdr3_fld)

    for seq in query:
        cdr3s += [seq.cdr3_len] * int(seq.copy_number)
    if len(cdr3s) == 0:
        return None, None
    q25, q75 = np.percentile(cdr3s, [25, 75])
    iqr = q75 - q25
    return float(q25 - 1.5 * iqr), float(q75 + 1.5 * iqr)


def _queue_tasks(session, sample_id, force, tasks):
    logger.info('Creating task queue to generate stats for sample {}.'.format(
        sample_id))
    existing_seq = session.query(Sequence).filter(
        Sequence.sample_id == sample_id)
    existing_nores = session.query(NoResult).filter(
        NoResult.sample_id == sample_id)
    if existing_seq.first() is None and existing_nores.first() is None:
        logger.warning('\tSKIPPING since there are no sequences in the '
                       'sample')
        return

    existing = session.query(SampleStats.sample_id).filter(
        SampleStats.sample_id == sample_id).first() is not None
    if force and existing:
        logger.warning('\tFORCING regeneration of stats')
    elif not force and existing:
        logger.warning('\tSKIPPING stats since they already exists and the '
                       '--force flag was not specified.')
        return

    min_cdr3, max_cdr3 = _get_cdr3_bounds(session, sample_id)
    for include_outliers in [True, False]:
        for only_full_reads in [True, False]:
            tasks.add_task({
                'func': 'seq',
                'sample_id': sample_id,
                'min_cdr3': min_cdr3,
                'max_cdr3': max_cdr3,
                'include_outliers': include_outliers,
                'only_full_reads': only_full_reads
            })
            tasks.add_task({
                'func': 'clone',
                'sample_id': sample_id,
                'min_cdr3': min_cdr3,
                'max_cdr3': max_cdr3,
                'include_outliers': include_outliers,
                'only_full_reads': only_full_reads
            })


def run_sample_stats(session, args):
    np.seterr(all='raise')
    mod_log.make_mod('sample_stats', session=session, commit=True,
                     info=vars(args))

    if args.sample_ids is None:
        samples = [s.id for s in session.query(Sample.id)]
    else:
        samples = args.sample_ids

    if args.force:
        q = session.query(SampleStats).filter(
            SampleStats.sample_id.in_(samples))
        q.delete(synchronize_session=False)
        session.commit()

    tasks = concurrent.TaskQueue()
    for sample_id in samples:
        _queue_tasks(session, sample_id, args.force, tasks)

    for i in range(0, args.nproc):
        session = config.init_db(args.db_config)
        tasks.add_worker(SampleStatsWorker(session))

    tasks.start()
    session.commit()
    session.close()
