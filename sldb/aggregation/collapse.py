from sqlalchemy import and_, distinct
from sqlalchemy.sql import desc, exists, text

import dnautils
import sldb.common.config as config
from sldb.common.models import (Clone, Sample, Sequence, SequenceCollapse,
                                Subject)
import sldb.common.modification_log as mod_log
import sldb.util.concurrent as concurrent

class CollapseWorker(concurrent.Worker):
    """A worker for collapsing sequences without including positions where
    either sequences has an 'N'.
    :param Session session: The database session
    """
    def __init__(self, session):
        self._session = session
        self._tasks = 0

    def do_task(self, bucket_hash):
        seqs = self._session.query(
            Sequence.sample_id, Sequence.ai, Sequence.sequence,
            Sequence.copy_number
        ).filter(
            Sequence.bucket_hash == bucket_hash
        ).all()

        to_process = sorted([{
            'sample_id': s.sample_id,
            'ai': s.ai,
            'sequence': s.sequence,
            'cn': s.copy_number
        } for s in seqs], key=lambda e: -e['cn'])

        while len(to_process) > 0:
            # Get the largest sequence in the list
            larger = to_process.pop(0)
            # Iterate over all smaller sequences to find matches
            for i in reversed(range(0, len(to_process))):
                smaller = to_process[i]
                if dnautils.equal(larger['sequence'], smaller['sequence']):
                    # Add the smaller sequence's copy number to the larger
                    larger['cn'] += smaller['cn']
                    # If the smaller sequence matches the larger, collapse it
                    # to the larger
                    self._session.add(SequenceCollapse(**{
                        'sample_id': smaller['sample_id'],
                        'seq_ai': smaller['ai'],
                        'copy_number_in_subject': 0,
                        'collapse_to_subject_seq_ai': larger['ai']
                    }))
                    # Delete the smaller sequence from the list to process
                    # since it's been collapsed
                    del to_process[i]

            # Update the larger sequence's copy number and "collapse" to itself
            self._session.add(SequenceCollapse(**{
                'sample_id': larger['sample_id'],
                'seq_ai': larger['ai'],
                'copy_number_in_subject': larger['cn'],
                'collapse_to_subject_seq_ai': larger['ai']
            }))

        self._tasks += 1
        if self._tasks % 100 == 0:
            self._session.commit()
            self._print('Collapsed {} buckets'.format(self._tasks))

    def cleanup(self):
        self._print('Committing collapsed sequences')
        self._session.commit()
        self._session.close()

def run_collapse(session, args):
    subject_ids = args.subject_ids or map(
        lambda e: e.id, session.query(Subject.id).all()
    )

    to_set_status = []
    for subject in subject_ids:
        if session.query(Sample).filter(
                Sample.subject_id == subject,
                Sample.status == 'identified').first() is None:
            print 'Subject {} already collapsed.  Skipping.'.format(subject)
            subject_ids.remove(subject)
        else:
            print 'Resetting collapse info for subject {}'.format(subject)
            samples = session.query(Sample).filter(
                  Sample.subject_id == subject).all()
            to_set_status.extend(samples)
            session.query(SequenceCollapse).filter(
                SequenceCollapse.sample_id.in_(map(lambda e: e.id, samples))
            ).delete(synchronize_session=False)
    session.commit()

    print 'Creating task queue to collapse {} subjects.'.format(
        len(subject_ids))

    tasks = concurrent.TaskQueue()

    for subject_id in subject_ids:
        buckets = session.query(
            Sequence.bucket_hash
        ).filter(
            Sequence.subject_id == subject_id
        ).group_by(
            Sequence.bucket_hash
        )
        for bucket in buckets:
            tasks.add_task(bucket.bucket_hash)

    for i in range(0, args.nproc):
        tasks.add_worker(CollapseWorker(config.init_db(args.db_config)))
    tasks.start()

    for sample in to_set_status:
        sample.status = 'collapsed'
    session.commit()
    session.close()
