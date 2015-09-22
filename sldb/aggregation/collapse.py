from sqlalchemy import and_, distinct
from sqlalchemy.sql import desc, exists, text

import dnautils
import sldb.common.config as config
from sldb.common.models import Clone, Sample, Sequence, Subject
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

    def do_task(self, args):
        """Initiates the correct method based on the type of collapsing that is
        being done
        :param dict args: A dictionary with at least a ``type`` key
        """
        if args['type'] == 'sample':
            self._collapse_sample(args)
        elif args['type'] == 'subject':
            self._collapse_subject(args)

    def _collapse_sample(self, args):
        """Collapses sequences in sample ``args['sample_id']`` with bucket hash
        ``args['bucket_hash``
        :param dict args: A dictionary with the keys ``sample_id``,
            ``bucket_hash``
        """
        seqs = self._session.query(
            Sequence.pk, Sequence.sequence, Sequence.copy_number
        ).filter(
            Sequence.sample_id == args['sample_id'],
            Sequence.bucket_hash == args['bucket_hash'],
        ).all()

        collapse_seqs(
            self._session, seqs, 'copy_number', 'copy_number_in_sample',
            'collapse_to_sample_seq_pk'
        )
        self._tasks += 1
        if self._tasks % 100 == 0:
            self._session.commit()
            self._print('Collapsed {} buckets'.format(self._tasks))

    def _collapse_subject(self, args):
        """Collapses all clones in the subject with ID ``args['subject_id']``
        and bucket ``(v_gene, j_gene, cdr3_num_nts)``
        :param dict args: A dictionary with the keys ``subject_id``,
            ``v_gene``, ``j_gene``, and ``cdr3_num_nts``
        """

        seqs = self._session.query(
            Sequence.pk, Sequence.sequence, Sequence.copy_number_in_sample
        ).filter(
            Sequence.subject_id == args['subject_id'],
            Sequence.bucket_hash == args['bucket_hash'],
            Sequence.copy_number_in_sample > 0,
        ).all()

        collapse_seqs(
            self._session, seqs, 'copy_number_in_sample',
            'copy_number_in_subject', 'collapse_to_subject_seq_pk',
        )

        self._tasks += 1
        if self._tasks % 100 == 0:
            self._session.commit()
            self._print('Collapsed {} buckets'.format(self._tasks))

    def cleanup(self):
        self._print('Committing collapsed sequences')
        self._session.commit()
        self._session.close()


def collapse_seqs(session, seqs, copy_field, collapse_copy_field,
                  collapse_pk_field):
    """Collapses sequences, aggregating their copy number in ``copy_field`` into
    a single sequences ``collapse_copy_field`` and setting the collapsed
    sequences' ``collapse_seq_id_field`` and ``collapse_sample_id_field`` to
    that of the collapsed-to sequence.
    :param Session session: The database session
    :param Sequence seqs: A list of sequences to collapse in any order
        containing at least the fields ``sample_id``, ``seq_id``, ``sequence``,
        and the fields with names passed in ``collapse_copy_field``,
        ``collapse_seq_id_field``, ``collapse_sample_id_field``.
    :param str copy_field: The name of the field which will be collapsed.  This
        field will be set to 0 for all sequences which are collapsed to
        another.
    :param str collapse_copy_field: The field to which ``copy_field`` values
        will be collapsed
    :param str collapse_seq_pk_field: The field to set to the ``pk`` of the
        sequence to which a sequence is collapsed
    """
    to_process = sorted([{
        'pk': s.pk,
        'sequence': s.sequence,
        'cn': getattr(s, copy_field)
    } for s in seqs], key=lambda e: -e['cn'])

    while len(to_process) > 0:
        # Get the largest sequence in the list
        larger = to_process[0]
        # Remove it from the list to process
        to_process = to_process[1:]
        # Iterate over all smaller sequences to find matches
        for i in reversed(range(0, len(to_process))):
            smaller = to_process[i]
            if dnautils.equal(larger['sequence'], smaller['sequence']):
                # Add the smaller sequence's copy number to the larger
                larger['cn'] += smaller['cn']
                # If the smaller sequence matches the larger, collapse it to
                # the larger
                session.query(Sequence).filter(
                    Sequence.pk == smaller['pk'],
                ).update({
                    collapse_pk_field: larger['pk'],
                    collapse_copy_field: 0
                })
                # Delete the smaller sequence from the list to process since
                # it's been collapsed
                del to_process[i]

        # Update the larger sequence's copy number and "collapse" to itself
        session.query(Sequence).filter(
            Sequence.pk == larger['pk'],
        ).update({
            collapse_pk_field: larger['pk'],
            collapse_copy_field: larger['cn']
        })


def collapse_samples(db_config, sample_ids, nproc):
    session = config.init_db(db_config)
    tasks = concurrent.TaskQueue()

    buckets = session.query(
        Sequence.bucket_hash, Sequence.sample_id
    ).group_by(
        Sequence.bucket_hash, Sequence.sample_id
    )
    for bucket in buckets:
        tasks.add_task({
            'type': 'sample',
            'sample_id': bucket.sample_id,
            'bucket_hash': bucket.bucket_hash
        })
    session.close()

    for i in range(0, nproc):
        worker_session = config.init_db(db_config)
        tasks.add_worker(CollapseWorker(worker_session))
    tasks.start()


def collapse_subjects(db_config, subject_ids, nproc):
    session = config.init_db(db_config)
    print 'Creating task queue to collapse {} subjects.'.format(
        len(subject_ids))

    tasks = concurrent.TaskQueue()

    for subject_id in subject_ids:
        buckets = session.query(
            Sequence.bucket_hash
        ).filter(
            Sequence.sample.has(subject_id=subject_id),
        ).group_by(
            Sequence.bucket_hash
        )
        for bucket in buckets:
            tasks.add_task({
                'type': 'subject',
                'subject_id': subject_id,
                'bucket_hash': bucket.bucket_hash
            })

    session.close()
    for i in range(0, nproc):
        worker_session = config.init_db(db_config)
        tasks.add_worker(CollapseWorker(worker_session))
    tasks.start()


def run_collapse(session, args):
    mod_log.make_mod('collapse', session=session, commit=True,
                     info=vars(args))

    subjects = args.subject_ids or map(
        lambda e: e.id, session.query(Subject.id).all()
    )

    potential_samples = session.query(
        Sample,
        exists().where(
            and_(
                Sequence.sample_id == Sample.id,
                Sequence.copy_number_in_sample > 0
            )
        )
    ).filter(
        Sample.subject_id.in_(subjects)
    ).all()

    if len(filter(lambda e: not e[1], potential_samples)) == 0:
        print 'All samples for subjects {} are already collapsed'.format(
            ','.join(map(str, subjects))
        )
        return

    to_collapse = []
    print 'Samples to process:'
    for sample, already_collapsed in potential_samples:
        if already_collapsed:
            print ('\tSample {} is not new.  Erasing old subject collapse '
                   'information.  Will re-collapse.').format(sample.id)
            session.query(Sequence).filter(
                Sequence.sample_id == sample.id
            ).update({
                'collapse_to_subject_seq_id': None,
                'collapse_to_subject_sample_id': None,
                'copy_number_in_subject': 0
            })
        else:
            print ('\tSample {} is new.  Collapsing at sample and subject '
                   'level.').format(sample.id)
            to_collapse.append(sample.id)

    session.close()
    print 'Collapsing {} samples'.format(len(to_collapse))
    collapse_samples(args.db_config, to_collapse, args.nproc)
    print 'Collapsing {} subjects'.format(len(subjects))
    collapse_subjects(args.db_config, subjects, args.nproc)
