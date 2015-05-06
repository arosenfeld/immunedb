import dnautils

from sqlalchemy import and_, distinct
from sqlalchemy.sql import desc, exists, text

import sldb.common.config as config
from sldb.common.models import Clone, CloneGroup, Sample, Sequence, Subject
import sldb.common.modification_log as mod_log
import sldb.util.concurrent as concurrent


class CollapseWorker(concurrent.Worker):
    def __init__(self, session):
        self._session = session
        self._tasks = 0

    def do_task(self, worker_id, args):
        if args['type'] == 'sample':
            self._collapse_sample(worker_id, args)
        elif args['type'] == 'subject':
            self._collapse_subject(worker_id, args)
        elif args['type'] == 'clone':
            self._collapse_clone(worker_id, args)

    def _collapse_sample(self, worker_id, args):
        seqs = self._session.query(
            Sequence.sample_id, Sequence.seq_id, Sequence.sequence,
            Sequence.copy_number
        ).filter(
            Sequence.sample_id == args['sample_id'],
            Sequence.v_gene == args['v_gene'],
            Sequence.j_gene == args['j_gene'],
            Sequence.junction_num_nts == args['junction_num_nts']
        ).all()
        self._print(worker_id, 'Collapsing bucket in sample {} with {} '
                    'seqs'.format(args['sample_id'], len(seqs)))

        collapse_seqs(
            self._session, seqs, 'copy_number', 'copy_number_in_sample',
            'collapse_to_sample_seq_id'
        )
        self._tasks += 1
        if self._tasks % 100 == 0:
            self._session.commit()

    def _collapse_subject(self, worker_id, args):
        seqs = self._session.query(
            Sequence.sample_id, Sequence.seq_id, Sequence.sequence,
            Sequence.copy_number_in_sample
        ).filter(
            Sequence.sample.has(subject_id=args['subject_id']),
            Sequence.copy_number_in_sample > 0,

            Sequence.v_gene == args['v_gene'],
            Sequence.j_gene == args['j_gene'],
            Sequence.junction_num_nts == args['junction_num_nts']
        ).all()
        self._print(worker_id, 'Collapsing bucket in subject {} with {} '
                    'seqs'.format(args['subject_id'], len(seqs)))

        collapse_seqs(
            self._session, seqs, 'copy_number_in_sample',
            'copy_number_in_subject', 'collapse_to_subject_seq_id',
            'collapse_to_subject_sample_id'
        )

        self._tasks += 1
        if self._tasks % 100 == 0:
            self._session.commit()

    def _collapse_clone(self, worker_id, args):
        self._print(worker_id, 'Collapsing clones {}'.format(
            ','.join(map(str, args['cids']))))

        for cid in args['cids']:
            seqs = self._session.query(
                Sequence.sample_id,
                Sequence.seq_id,
                Sequence.sequence,
                Sequence.copy_number_in_subject
            ).filter(
                Sequence.clone_id == cid,
                Sequence.copy_number_in_subject > 0
            ).all()

            collapse_seqs(
                self._session, seqs, 'copy_number_in_subject', 'copy_number_in_clone',
                'collapse_to_clone_seq_id', 'collapse_to_clone_sample_id'
            )

        self._tasks += 1
        if self._tasks % 100 == 0:
            self._session.commit()


    def cleanup(self, worker_id):
        self._print(worker_id, 'Committing collapsed sequences')
        self._session.commit()
        self._session.close()


def collapse_seqs(session, seqs, copy_field, collapse_copy_field,
                  collapse_seq_id_field, collapse_sample_id_field=None):
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
    print 'Creating task queue to collapse {} samples.'.format(len(sample_ids))

    session.query(Sequence).filter(
        Sequence.sample_id.in_(sample_ids)
    ).update({
        'collapse_to_subject_seq_id': None,
        'collapse_to_subject_sample_id': None,
        'copy_number_in_subject': 0,
        'collapse_to_clone_seq_id': None,
        'collapse_to_clone_sample_id': None,
        'copy_number_in_clone': 0,
    }, synchronize_session=False)
    session.commit()

    tasks = concurrent.TaskQueue()

    for sample_id in sample_ids:
        buckets = session.query(
            Sequence.v_gene, Sequence.j_gene, Sequence.junction_num_nts,
        ).filter(
            Sequence.sample_id == sample_id
        ).group_by(
            Sequence.v_gene, Sequence.j_gene, Sequence.junction_num_nts,
        )
        for bucket in buckets:
            tasks.add_task({
                'type': 'sample',
                'sample_id': sample_id,
                'v_gene': bucket.v_gene,
                'j_gene': bucket.j_gene,
                'junction_num_nts': bucket.junction_num_nts
            })
    session.close()

    for i in range(0, nproc):
        worker_session = config.init_db(master_db_config, data_db_config)
        tasks.add_worker(CollapseWorker(worker_session))
    tasks.start()


def collapse_subjects(master_db_config, data_db_config, subject_ids,
                      nproc):
    session = config.init_db(master_db_config, data_db_config)
    print 'Creating task queue to collapse {} subjects.'.format(
        len(subject_ids))

    session.query(Sequence).filter(
        exists().where(
            and_(
                Sample.id == Sequence.sample_id,
                Sample.subject_id.in_(subject_ids)
            )
        )
    ).update({
        'collapse_to_clone_seq_id': None,
        'collapse_to_clone_sample_id': None,
        'copy_number_in_clone': 0,
    }, synchronize_session=False)
    session.commit()

    tasks = concurrent.TaskQueue()

    for subject_id in subject_ids:
        buckets = session.query(
            Sequence.v_gene, Sequence.j_gene, Sequence.junction_num_nts
        ).filter(
            Sequence.sample.has(subject_id=subject_id),
        ).group_by(
            Sequence.v_gene, Sequence.j_gene, Sequence.junction_num_nts
        )
        for bucket in buckets:
            tasks.add_task({
                'type': 'subject',
                'subject_id': subject_id,
                'v_gene': bucket.v_gene,
                'j_gene': bucket.j_gene,
                'junction_num_nts': bucket.junction_num_nts
            })

    session.close()
    for i in range(0, nproc):
        worker_session = config.init_db(master_db_config,
                                        data_db_config)
        tasks.add_worker(CollapseWorker(worker_session))
    tasks.start()


def collapse_clones(master_db_config, data_db_config, clone_ids,
                    nproc, chunk_size=25):
    '''
    print 'Creating task queue to collapse {} clones.'.format(len(clone_ids))
    tasks = concurrent.TaskQueue()
    # Collapsing clones is relatively fast and the overhead of context switching
    # is high, so give each worker a few clones to collapse instead of just one.
    for i in range(0, len(clone_ids), chunk_size):
        tasks.add_task({
            'type': 'clone',
            'cids': clone_ids[i:i+chunk_size]
        })

    for i in range(0, nproc):
        worker_session = config.init_db(master_db_config, data_db_config)
        tasks.add_worker(CollapseWorker(worker_session))

    tasks.start()
    '''

    session = config.init_db(master_db_config, data_db_config)
    print 'Pushing clone IDs to subject sequences'
    session.connection(mapper=Sequence).execute(text('''
        update
            sequences as s
        join
            (select seq_id, sample_id, clone_id from sequences where
                copy_number_in_clone > 0) as j
        on
            (s.collapse_to_clone_seq_id=j.seq_id and
                s.collapse_to_clone_sample_id=j.sample_id)
        set
            s.clone_id=j.clone_id
        where
            s.copy_number_in_clone=0
    '''))

    print 'Pushing clone IDs to sample sequences'
    session.connection(mapper=Sequence).execute(text('''
        update
            sequences as s
        join
            (select seq_id, sample_id, clone_id from sequences where
                copy_number_in_subject > 1) as j
        on
            (s.collapse_to_subject_seq_id=j.seq_id and
                s.collapse_to_subject_sample_id=j.sample_id)
        set
            s.clone_id=j.clone_id
        where
            s.copy_number_in_clone=0
    '''))

    print 'Pushing clone IDs to individual sequences'
    session.connection(mapper=Sequence).execute(text('''
        update
            sequences as s
        join
            (select seq_id, sample_id, clone_id from sequences where
                copy_number_in_sample > 0) as j
        on
            (s.collapse_to_sample_seq_id=j.seq_id and
                s.sample_id=j.sample_id)
        set
            s.clone_id=j.clone_id
        where
            s.copy_number_in_sample=0
    '''))
    session.commit()
    session.close()


def run_collapse(session, args):
    mod_log.make_mod('collapse', session=session, commit=True,
                     info=vars(args))
    if args.level == 'samples':
        ids = args.ids or map(
            lambda r: r.id, session.query(Sample).all()
        )
        func = collapse_samples
    elif args.level == 'subjects':
        ids = args.ids or map(
            lambda r: r.id, session.query(Subject).all()
        )
        func = collapse_subjects
    elif args.level == 'clones':
        func = collapse_clones
        if args.ids is not None:
            ids = map(lambda r: r.id, session.query(Clone.id).filter(
                Clone.subject_id.in_(args.ids)).all()
            )
        else:
            ids = map(lambda r: r.id, session.query(Clone).all())
    else:
        print '[ERROR] Unknown collapse level.'
        return
    session.close()
    func(args.master_db_config, args.data_db_config,
         ids, args.nproc)
