import json

from sqlalchemy import distinct, func
from sqlalchemy.orm import scoped_session

import sldb.common.baseline as baseline
import sldb.common.config as config
from sldb.common.models import Clone, CloneStats, Sequence
import sldb.common.modification_log as mod_log
from sldb.common.mutations import CloneMutations
import sldb.util.concurrent as concurrent


class CloneStatsWorker(concurrent.Worker):
    def __init__(self, session, baseline_path, baseline_temp):
        self._session = session
        self._baseline_path = baseline_path
        self._baseline_temp = baseline_temp
        self._total_completed = 0

    def do_task(self, worker_id, args):
        clone_id = args['clone_id']
        sample_id = args['sample_id']
        if sample_id != 0:
            self._sample_stats(worker_id, clone_id, sample_id)
        else:
            self._total_stats(worker_id, clone_id)

    def _sample_stats(self, worker_id, clone_id, sample_id):
        self._print(worker_id, 'Clone {}, sample {}'.format(
            clone_id, sample_id))

        existing = self._session.query(CloneStats).filter(
            CloneStats.clone_id == clone_id,
            CloneStats.sample_id == sample_id).first()
        if existing is not None:
            return

        counts = self._session.query(
            func.count(distinct(Sequence.sequence_replaced)).label('unique'),
            func.sum(Sequence.copy_number).label('total')
        ).filter(
            Sequence.clone_id == clone_id,
            Sequence.sample_id == sample_id
        ).first()

        sample_mutations = CloneMutations(
            self._session,
            self._session.query(Clone).filter(
                Clone.id == clone_id).first()
        ).calculate(
            commit_seqs=True, limit_samples=[sample_id],
            mode=CloneMutations.MODE_SAMPLES_ONLY
        )[sample_id]

        selection_pressure = baseline.get_selection(
            self._session, clone_id, self._baseline_path,
            samples=[sample_id], temp_dir=self._baseline_temp)

        self._session.add(CloneStats(
            clone_id=clone_id,
            sample_id=sample_id,
            unique_cnt=counts.unique,
            total_cnt=counts.total,
            mutations=json.dumps(sample_mutations.get_all()),
            selection_pressure=json.dumps(selection_pressure)
        ))

        self._total_completed += 1
        if self._total_completed % 1000 == 0:
            self._session.commit()

    def _total_stats(self, worker_id, clone_id):
        existing = self._session.query(CloneStats).filter(
            CloneStats.clone_id == clone_id,
            CloneStats.sample_id == 0).first()
        if existing is not None:
            return

        self._print(worker_id, 'Clone {}, all samples'.format(clone_id))
        # Get the counts for the entire clone
        counts = self._session.query(
            func.count(distinct(Sequence.sequence_replaced)).label('unique'),
            func.sum(Sequence.copy_number).label('total')
        ).filter(Sequence.clone_id == clone_id).first()

        total_mutations = CloneMutations(
            self._session,
            self._session.query(Clone).filter(
                Clone.id == clone_id).first()
        ).calculate(mode=CloneMutations.MODE_TOTAL_ONLY)

        # Get selection pressure with Baseline
        selection_pressure = baseline.get_selection(
            self._session, clone_id, self._baseline_path,
            temp_dir=self._baseline_temp)

        # Add the statistics for the whole clone, denoted with a 0 in
        # the sample_id field
        self._session.add(CloneStats(
            clone_id=clone_id,
            sample_id=0,
            unique_cnt=counts.unique,
            total_cnt=counts.total,
            mutations=json.dumps(total_mutations.get_all()),
            selection_pressure=json.dumps(selection_pressure)
        ))

    def cleanup(self, worker_id):
        self._session.commit()


def run_clone_stats(session, args):
    mod_log.make_mod('clone_stats', session=session, commit=True,
                     info=vars(args))

    if args.clone_ids is not None:
        clones = args.clone_ids
    elif args.subjects is not None:
        clones = map(lambda c: c.id, session.query(Clone.id).filter(
            Clone.subject_id.in_(args.subjects)).all())
    else:
        clones = map(lambda c: c.id, session.query(Clone.id).all())

    if args.force:
        session.query(CloneStats).filter(
            CloneStats.clone_id.in_(clones)
        ).delete(synchronize_session=False)
        session.commit()

    tasks = concurrent.TaskQueue()
    for cid in clones:
        tasks.add_task({
            'clone_id': cid,
            'sample_id': 0
        })
        for sid in map(lambda c: c.sample_id, session.query(
                distinct(Sequence.sample_id).label('sample_id')
                ).filter(Sequence.clone_id == cid)):
            tasks.add_task({
                'clone_id': cid,
                'sample_id': sid
            })

    for i in range(0, args.nproc):
        session = config.init_db(args.master_db_config, args.data_db_config)
        tasks.add_worker(CloneStatsWorker(
            session, args.baseline_path, args.temp))

    tasks.start()
