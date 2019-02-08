import json

from sqlalchemy import distinct, func

import immunedb.common.config as config
from immunedb.common.models import (Clone, CloneStats, Sequence,
                                    SequenceCollapse)
import immunedb.common.modification_log as mod_log
from immunedb.common.mutations import CloneMutations
import immunedb.util.concurrent as concurrent
from immunedb.util.log import logger


class CloneStatsWorker(concurrent.Worker):
    """A worker class for generating clone statistics.  This worker will accept
    one clone at a time for parallelization.

    :param Session session: The database session

    """
    def __init__(self, session):
        self._session = session

    def do_task(self, clone_id):
        """Starts the task of generating clone statistics for a given
        clone_id.
        :param int args: The clone_id for which to calculate statistics

        """

        existing = self._session.query(CloneStats).filter(
            CloneStats.clone_id == clone_id).first()

        if existing is not None:
            return

        self.info('Clone {}'.format(clone_id))
        sample_ids = [c.sample_id for c in self._session.query(
                distinct(Sequence.sample_id).label('sample_id')
            ).filter(
                Sequence.clone_id == clone_id
            )]
        sample_ids.append(None)
        for sample_id in sample_ids:
            self._process_sample(clone_id, sample_id)

    def _process_sample(self, clone_id, sample_id):
        """Processes clone statistics for one sample (or the aggregate of all
        samples).  If ``sample_id`` is None the statistics for all sequences in
        the clone is generated.

        :param int clone_id: The ID of the clone
        :param int sample_id: The ID of a sample in which the clone exists

        """

        top_seq = self._session.query(
            Sequence.ai,
            Sequence.copy_number,
            Sequence.sequence
        ).filter(
            Sequence.clone_id == clone_id
        )
        if sample_id is None:
            counts = self._session.query(
                func.count(Sequence.ai).label('unique'),
                func.sum(SequenceCollapse.copy_number_in_subject).label(
                    'total'),
                func.sum(
                    (Sequence.v_match / Sequence.v_length) *
                    Sequence.copy_number
                ).label('v_identity')
            ).join(SequenceCollapse).filter(
                Sequence.clone_id == clone_id,
                SequenceCollapse.copy_number_in_subject > 0
            ).first()
        else:
            counts = self._session.query(
                func.count(Sequence.ai).label('unique'),
                func.sum(Sequence.copy_number).label('total'),
                func.sum(
                    (Sequence.v_match / Sequence.v_length) *
                    Sequence.copy_number
                ).label('v_identity')
            ).filter(
                Sequence.sample_id == sample_id,
                Sequence.clone_id == clone_id
            ).first()
            top_seq = top_seq.filter(Sequence.sample_id == sample_id)

        top_seq = top_seq.order_by(Sequence.copy_number.desc()).first()

        clone_inst = self._session.query(Clone).filter(
            Clone.id == clone_id).first()
        sample_mutations = CloneMutations(
            self._session,
            clone_inst
        ).calculate(
            commit_seqs=sample_id is not None, limit_samples=[sample_id],
        )[sample_id]

        record_values = {
            'clone_id': clone_id,
            'subject_id': clone_inst.subject_id,
            'functional': clone_inst.functional,
            'unique_cnt': counts.unique,
            'total_cnt': counts.total,
            'mutations': json.dumps(sample_mutations.get_all()),
            'avg_v_identity': counts.v_identity / counts.total,
            'top_copy_seq_ai': top_seq.ai,
            'top_copy_seq_sequence': top_seq.sequence,
            'top_copy_seq_copies': top_seq.copy_number
        }

        self._session.add(CloneStats(sample_id=sample_id, **record_values))

        # If this clone only appears in one sample, the 'total clone' stats are
        # the same as for the single sample
        if sample_id is None:
            clone_inst.overall_unique_cnt = counts.unique
            clone_inst.overall_instance_cnt = self._session.query(
                func.count(Sequence.ai).label('instances'),
            ).filter(
                Sequence.clone_id == clone_id
            ).first().instances

            clone_inst.overall_total_cnt = counts.total
        self._session.commit()

    def cleanup(self):
        self._session.commit()
        self._session.close()


def run_clone_stats(session, args):
    """Runs the clone statistics generation stage of the pipeline.
    :param Session session: The database session
    :param Namespace args: The arguments passed to the command

    """
    mod_log.make_mod('clone_stats', session=session, commit=True,
                     info=vars(args))

    if args.clone_ids is not None:
        clones = args.clone_ids
    elif args.subject_ids is not None:
        clones = [c.id for c in session.query(Clone.id).filter(
            Clone.subject_id.in_(args.subject_ids))]
    else:
        clones = [c.id for c in session.query(Clone.id)]
    clones.sort()

    if args.regen:
        logger.info('Deleting old clone statistics for {} clones'.format(
            len(clones)))
        session.query(CloneStats).filter(
            CloneStats.clone_id.in_(clones)
        ).delete(synchronize_session=False)
        session.commit()

    tasks = concurrent.TaskQueue()
    logger.info('Creating task queue to generate stats for {} clones.'.format(
        len(clones)))
    for cid in clones:
        tasks.add_task(cid)

    for i in range(0, args.nproc):
        session = config.init_db(args.db_config)
        tasks.add_worker(CloneStatsWorker(session))

    tasks.start()
