from sqlalchemy import and_, distinct, or_
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
        """Collapses sequences in sample.

        :param dict args: A dictionary with the key ``args['sample_id']``

        """
        buckets = self._session.query(
            Sequence.v_gene, Sequence.j_gene, Sequence.cdr3_num_nts,
        ).filter(
            Sequence.sample_id == args['sample_id']
        ).group_by(
            Sequence.v_gene, Sequence.j_gene, Sequence.cdr3_num_nts,
        )
        num_buckets = buckets.count()

        for i, bucket in enumerate(buckets):
            seqs = self._session.query(
                Sequence.sample_id, Sequence.seq_id, Sequence.sequence,
                Sequence.copy_number
            ).filter(
                Sequence.sample_id == args['sample_id'],
                Sequence.v_gene == bucket.v_gene,
                Sequence.j_gene == bucket.j_gene,
                Sequence.cdr3_num_nts == bucket.cdr3_num_nts
            ).all()

            collapse_seqs(
                self._session, seqs, 'copy_number', 'copy_number_in_sample',
                'collapse_to_sample_seq_id'
            )
            if i % 100 == 0:
                self._print('\tCollapsed {}/{} buckets in sample {}'.format(
                    i, num_buckets, args['sample_id']))
                self._session.commit()
        self._session.commit()
        self._print('\tFinished collapsing sample {}'.format(
            args['sample_id']))

    def _collapse_subject(self, args):
        """Collapses all clones in the subject with ID ``args['subject_id']``

        :param dict args: A dictionary with the key ``subject_id``

        """

        buckets = self._session.query(
            Sequence.v_gene, Sequence.j_gene, Sequence.cdr3_num_nts
        ).filter(
            Sequence.sample.has(subject_id=args['subject_id']),
        ).group_by(
            Sequence.v_gene, Sequence.j_gene, Sequence.cdr3_num_nts
        )
        num_buckets = buckets.count()

        for i, bucket in enumerate(buckets):
            seqs = self._session.query(
                Sequence.sample_id, Sequence.seq_id, Sequence.sequence,
                Sequence.copy_number_in_sample
            ).filter(
                Sequence.sample.has(subject_id=args['subject_id']),
                Sequence.copy_number_in_sample > 0,

                Sequence.v_gene == bucket.v_gene,
                Sequence.j_gene == bucket.j_gene,
                Sequence.cdr3_num_nts == bucket.cdr3_num_nts
            ).all()

            collapse_seqs(
                self._session, seqs, 'copy_number_in_sample',
                'copy_number_in_subject', 'collapse_to_subject_seq_id',
                'collapse_to_subject_sample_id'
            )

            if i % 100 == 0:
                self._print('\tCollapsed {}/{} buckets in subject {}'.format(
                    i, num_buckets, args['subject_id']))
                self._session.commit()
        self._session.commit()
        self._print('\tFinished collapsing subject {}'.format(
            args['subject_id']))

    def cleanup(self):
        self._session.close()


def collapse_seqs(session, seqs, copy_field, collapse_copy_field,
                  collapse_seq_id_field, collapse_sample_id_field=None):
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
    :param str collapse_seq_id_field: The field to set to the ``seq_id`` of the
        sequence to which a sequence is collapsed
    :param str collapse_sample_id_field: The field to set to the ``sample_id``
        of the sequence to which a sequence is collapsed

    """
    to_process = sorted([{
        'sample_id': s.sample_id,
        'seq_id': s.seq_id,
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
                update_dict = {
                    collapse_seq_id_field: larger['seq_id'],
                    collapse_copy_field: 0
                }
                if collapse_sample_id_field is not None:
                    update_dict[collapse_sample_id_field] = larger['sample_id']
                session.query(Sequence).filter(
                    Sequence.sample_id == smaller['sample_id'],
                    Sequence.seq_id == smaller['seq_id']
                ).update(update_dict)
                # Delete the smaller sequence from the list to process since
                # it's been collapsed
                del to_process[i]

        # Update the larger sequence's copy number and "collapse" to itself
        update_dict = {
            collapse_seq_id_field: larger['seq_id'],
            collapse_copy_field: larger['cn']
        }
        if collapse_sample_id_field is not None:
            update_dict[collapse_sample_id_field] = larger['sample_id']
        session.query(Sequence).filter(
            Sequence.sample_id == larger['sample_id'],
            Sequence.seq_id == larger['seq_id'],
        ).update(update_dict)


def collapse_samples(master_db_config, data_db_config, sample_ids,
                     nproc):
    session = config.init_db(master_db_config, data_db_config)
    tasks = concurrent.TaskQueue()

    for sample_id in sample_ids:
        tasks.add_task({
            'type': 'sample',
            'sample_id': sample_id
        })
    session.close()

    for i in range(0, nproc):
        worker_session = config.init_db(master_db_config, data_db_config)
        tasks.add_worker(CollapseWorker(worker_session))
    tasks.start()


def collapse_subjects(master_db_config, data_db_config, subject_ids,
                      nproc):
    session = config.init_db(master_db_config, data_db_config)

    tasks = concurrent.TaskQueue()

    for subject_id in subject_ids:
        tasks.add_task({
            'type': 'subject',
            'subject_id': subject_id,
        })

    session.close()
    for i in range(0, nproc):
        worker_session = config.init_db(master_db_config,
                                        data_db_config)
        tasks.add_worker(CollapseWorker(worker_session))
    tasks.start()


def run_collapse(session, args):
    mod_log.make_mod('collapse', session=session, commit=True,
                     info=vars(args))

    subjects = args.subject_ids or map(
        lambda e: e.id, session.query(Subject.id).all()
    )

    samples = map(lambda e: e.id, session.query(
        Sample.id
    ).filter(
        Sample.subject_id.in_(subjects)).all()
    )
    for sample in samples:
        if session.query(Sequence).filter(
                Sequence.sample_id == sample,
                or_(Sequence.copy_number_in_sample > 0,
                    Sequence.copy_number_in_subject > 0)).first() is None:
            continue
        print 'Clearing old collapse information for sample {}'.format(sample)
        session.query(Sequence).filter(Sequence.sample_id == sample).filter(
        ).update({
            'collapse_to_sample_seq_id': None,
            'copy_number_in_sample': 0,

            'collapse_to_subject_seq_id': None,
            'collapse_to_subject_sample_id': None,
            'copy_number_in_subject': 0
        })
        session.commit()

    session.close()
    print 'Collapsing {} samples'.format(len(samples))
    collapse_samples(args.master_db_config, args.data_db_config, samples,
                     args.nproc)
    print 'Collapsing {} subjects'.format(len(subjects))
    collapse_subjects(args.master_db_config, args.data_db_config, subjects,
                     args.nproc)
