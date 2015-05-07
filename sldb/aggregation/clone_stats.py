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

    def do_task(self, worker_id, args):
        clone_id = args['clone_id']
        sample_id = args['sample_id']
        single = args['single']
        self._context_stats(worker_id, clone_id, sample_id, single)

    def _context_stats(self, worker_id, clone_id, sample_id, single):
        self._print(worker_id, 'Clone {}, sample {}'.format(
            clone_id, sample_id if sample_id != 0 else 'ALL')
        )

        existing = self._session.query(CloneStats).filter(
            CloneStats.clone_id == clone_id,
            CloneStats.sample_id == sample_id).first()

        if existing is not None:
            return

        if sample_id != 0:
            counts = self._session.query(
                func.count(Sequence.seq_id).label('unique'),
                func.sum(Sequence.copy_number_in_sample).label('total')
            ).filter(
                Sequence.clone_id == clone_id,
                Sequence.sample_id == sample_id,
                Sequence.copy_number_in_sample > 0
            ).first()

            sample_mutations = CloneMutations(
                self._session,
                self._session.query(Clone).filter(Clone.id == clone_id).first()
            ).calculate(
                commit_seqs=True, limit_samples=[sample_id],
            )[sample_id]
        else:
            counts = self._session.query(
                func.count(Sequence.seq_id).label('unique'),
                func.sum(Sequence.copy_number_in_clone).label('total')
            ).filter(
                Sequence.clone_id == clone_id,
                Sequence.copy_number_in_clone > 0
            ).first()

            sample_mutations = CloneMutations(
                self._session,
                self._session.query(Clone).filter(Clone.id == clone_id).first()
            ).calculate(limit_samples=[0])[0]

        selection_pressure = baseline.get_selection(
            self._session, clone_id, self._baseline_path,
            samples=[sample_id] if sample_id != 0 else None,
            temp_dir=self._baseline_temp)

        record_values = {
            'clone_id': clone_id,
            'unique_cnt': counts.unique,
            'total_cnt': counts.total,
            'mutations': json.dumps(sample_mutations.get_all()),
            'selection_pressure': json.dumps(selection_pressure)
        }
        self._session.add(CloneStats(sample_id=sample_id, **record_values))

        # If this clone only appears in one sample, the 'total clone' stats are
        # the same as for the single sample
        if single:
            self._session.add(CloneStats(sample_id=0, **record_values))


    def cleanup(self, worker_id):
        self._session.commit()
        self._session.close()


def run_clone_stats(session, args):
    mod_log.make_mod('clone_stats', session=session, commit=True,
                     info=vars(args))

    if args.clone_ids is not None:
        clones = args.clone_ids
    elif args.subject_ids is not None:
        clones = map(lambda c: c.id, session.query(Clone.id).filter(
            Clone.subject_id.in_(args.subject_ids)).all())
    else:
        clones = map(lambda c: c.id, session.query(Clone.id).all())

    if args.force:
        session.query(CloneStats).filter(
            CloneStats.clone_id.in_(clones)
        ).delete(synchronize_session=False)
        session.commit()

    tasks = concurrent.TaskQueue()
    print 'Creating task queue to generate stats for {} clones.'.format(
        len(clones)
    )
    for cid in clones:
        sample_ids = map(lambda c: c.sample_id, session.query(
            distinct(Sequence.sample_id).label('sample_id')
        ).filter(
            Sequence.clone_id == cid,
            Sequence.copy_number_in_sample > 0
        ))
        if len(sample_ids) > 1:
            tasks.add_task({
                'clone_id': cid,
                'single': False,
                'sample_id': 0
            })
            for sid in sample_ids:
                tasks.add_task({
                    'clone_id': cid,
                    'single': False,
                    'sample_id': sid
                })
        else:
            tasks.add_task({
                'clone_id': cid,
                'single': True,
                'sample_id': sample_ids[0]
            })


    for i in range(0, args.nproc):
        session = config.init_db(args.master_db_config, args.data_db_config)
        tasks.add_worker(CloneStatsWorker(
            session, args.baseline_path, args.temp))

    tasks.start()
