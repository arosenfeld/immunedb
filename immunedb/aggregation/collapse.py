from sqlalchemy.sql import exists

import dnautils
import immunedb.common.config as config
from immunedb.common.models import (Clone, Sample, Sequence, SequenceCollapse,
                                    Subject)
import immunedb.common.modification_log as mod_log
import immunedb.util.concurrent as concurrent

from immunedb.util.log import logger


class CollapseWorker(concurrent.Worker):
    """A worker for collapsing sequences without including positions where
    either sequences has an 'N'.
    :param Session session: The database session
    """
    def __init__(self, session):
        self._session = session
        self._tasks = 0

    def do_task(self, bucket):
        seqs = self._session.query(
            Sequence.sample_id, Sequence.ai, Sequence.seq_id,
            Sequence.sequence, Sequence.copy_number
        ).filter(
            Sequence.subject_id == bucket.subject_id,
            Sequence.v_gene == bucket.v_gene,
            Sequence.j_gene == bucket.j_gene,
            Sequence.cdr3_num_nts == bucket.cdr3_num_nts,
            Sequence._insertions == bucket._insertions,
            Sequence._deletions == bucket._deletions
        )

        to_process = sorted([{
            'sample_id': s.sample_id,
            'ai': s.ai,
            'seq_id': s.seq_id,
            'sequence': s.sequence,
            'cn': s.copy_number
        } for s in seqs], key=lambda e: -e['cn'])

        while len(to_process) > 0:
            # Get the largest sequence in the list
            larger = to_process.pop(0)
            # Iterate over all smaller sequences to find matches
            instances = 1
            samples = set([larger['sample_id']])
            for i in reversed(range(len(to_process))):
                smaller = to_process[i]
                if len(larger['sequence']) != len(smaller['sequence']):
                    self.warning('Tried to collapse sequences of different '
                                 'lengths.  AIs are {} {}'.format(
                                     larger['ai'], smaller['ai']))
                elif dnautils.equal(larger['sequence'], smaller['sequence']):
                    # Add the smaller sequence's copy number to the larger
                    larger['cn'] += smaller['cn']
                    # If the smaller sequence matches the larger, collapse it
                    # to the larger
                    self._session.add(SequenceCollapse(**{
                        'sample_id': smaller['sample_id'],
                        'seq_ai': smaller['ai'],
                        'collapse_to_subject_seq_ai': larger['ai'],
                        'collapse_to_subject_sample_id': larger['sample_id'],
                        'collapse_to_subject_seq_id': larger['seq_id'],
                        'instances_in_subject': 0,
                        'copy_number_in_subject': 0,
                        'samples_in_subject': 0,
                    }))
                    instances += 1
                    samples.add(smaller['sample_id'])
                    # Delete the smaller sequence from the list to process
                    # since it's been collapsed
                    del to_process[i]

            # Update the larger sequence's copy number and "collapse" to itself
            self._session.add(SequenceCollapse(**{
                'sample_id': larger['sample_id'],
                'seq_ai': larger['ai'],
                'collapse_to_subject_sample_id': larger['sample_id'],
                'collapse_to_subject_seq_id': larger['seq_id'],
                'collapse_to_subject_seq_ai': larger['ai'],
                'instances_in_subject': instances,
                'copy_number_in_subject': larger['cn'],
                'samples_in_subject': len(samples),
            }))

        self._session.commit()
        self._tasks += 1
        if self._tasks > 0 and self._tasks % 100 == 0:
            self.info('Collapsed {} buckets'.format(self._tasks))

    def cleanup(self):
        self.info('Committing collapsed sequences')
        self._session.commit()
        self._session.close()


def run_collapse(session, args):
    mod_log.make_mod('collapse', session=session, commit=True,
                     info=vars(args))
    subject_ids = []

    subjects = (args.subject_ids or [e.id for e in session.query(Subject.id)])
    for subject in subjects:
        if session.query(Sample).filter(
                Sample.subject_id == subject,
                ~exists().where(
                    SequenceCollapse.sample_id == Sample.id
                )).first() is None:
            logger.info('Subject {} already collapsed.  Skipping.'.format(
                subject))
        else:
            logger.info('Resetting collapse info for subject {}'.format(
                subject))
            samples = session.query(Sample).filter(
                  Sample.subject_id == subject
            )
            for sample in samples:
                session.query(SequenceCollapse).filter(
                    SequenceCollapse.sample_id == sample.id
                ).delete(synchronize_session=False)
                sample.sample_stats = []
            logger.info('Resetting clone info for subject {}'.format(subject))
            session.query(Clone).filter(Clone.subject_id == subject).delete()
            subject_ids.append(subject)
    session.commit()

    logger.info('Creating task queue to collapse {} subjects.'.format(
        len(subject_ids)))

    tasks = concurrent.TaskQueue()

    for subject_id in subject_ids:
        buckets = session.query(
            Sequence.subject_id, Sequence.v_gene, Sequence.j_gene,
            Sequence.cdr3_num_nts, Sequence._insertions, Sequence._deletions
        ).filter(
            Sequence.subject_id == subject_id
        ).group_by(
            Sequence.subject_id, Sequence.v_gene, Sequence.j_gene,
            Sequence.cdr3_num_nts, Sequence._insertions, Sequence._deletions
        )
        for bucket in buckets:
            tasks.add_task(bucket)

    logger.info('Generated {} total tasks'.format(tasks.num_tasks()))

    for i in range(0, min(tasks.num_tasks(), args.nproc)):
        tasks.add_worker(CollapseWorker(config.init_db(args.db_config)))
    tasks.start()

    session.close()
