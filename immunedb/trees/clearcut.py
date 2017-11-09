import json

import immunedb.common.config as config
from immunedb.common.models import Clone, Sequence, SequenceCollapse
import immunedb.common.modification_log as mod_log
from immunedb.trees import PhylogeneticTree, tree_as_dict
import immunedb.util.concurrent as concurrent
from immunedb.util.log import logger


class ClearcutWorker(concurrent.Worker):
    def __init__(self, session, tree_prog, min_count, min_samples,
                 min_seq_copies, exclude_stops):
        self._session = session
        self._tree_prog = tree_prog
        self._min_count = min_count
        self._min_samples = min_samples
        self._min_seq_copies = min_seq_copies
        self._exclude_stops = exclude_stops

    def do_task(self, clone_id):
        clone_inst = self._session.query(Clone).filter(
            Clone.id == clone_id).first()
        if clone_inst is None:
            return

        self.info('Running clone {}'.format(clone_inst.id))

        sequences = self._session.query(
            Sequence
        ).join(SequenceCollapse).filter(
            Sequence.clone_id == clone_id,
            SequenceCollapse.copy_number_in_subject > self._min_seq_copies
        )

        if self._exclude_stops:
            sequences = sequences.filter(Sequence.stop == 0)

        sequences = sequences.order_by(Sequence.v_length)

        try:
            tree = PhylogeneticTree(clone_inst.consensus_germline, sequences)
            tree.run(self._session, self._tree_prog)
        except Exception as e:
            logger.error('Error running clone {}: {}'.format(clone_id, e))
            return

        final = {
            'info': {
                'min_count': self._min_count,
                'min_samples': self._min_samples,
                'min_seq_copies': self._min_seq_copies,
                'exclude_stops': self._exclude_stops
            },
            'tree': tree_as_dict(tree.tree)
        }
        clone_inst.tree = json.dumps(final)
        self._session.add(clone_inst)
        self._session.commit()


def run_clearcut(session, args):
    if args.clone_ids is not None:
        clones = session.query(Clone.id).filter(
            Clone.id.in_(args.clone_ids))
    else:
        if args.subject_ids is not None:
            clones = session.query(Clone.id).filter(
                Clone.subject_id.in_(args.subject_ids))
        else:
            clones = session.query(Clone.id)

    if not args.force:
        clones = clones.filter(Clone.tree.is_(None))
    clones = [c.id for c in clones]
    mod_log.make_mod('clone_tree', session=session, commit=True,
                     info=vars(args))

    tasks = concurrent.TaskQueue()

    logger.info('Creating task queue for clones')
    for clone_id in clones:
        tasks.add_task(clone_id)

    for _ in range(0, args.nproc):
        session = config.init_db(args.db_config)
        tasks.add_worker(ClearcutWorker(session, args.clearcut_path,
                                        args.min_count, args.min_samples,
                                        args.min_seq_copies,
                                        args.exclude_stops))

    tasks.start()
