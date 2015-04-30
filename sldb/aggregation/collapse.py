import re

from sqlalchemy import and_, distinct
from sqlalchemy.sql import desc, exists, text

import sldb.common.config as config
from sldb.common.models import Clone, Sample, Sequence, Subject
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

    def _collapse_sample(self, worker_id, args):
        seqs = self._session.query(
            Sequence.sample_id, Sequence.seq_id, Sequence.sequence,
            Sequence.copy_number
        ).filter(
            Sequence.sample_id == args['sample_id'],
            Sequence.v_gene == args['v_gene'],
            Sequence.j_gene == args['j_gene'],
            Sequence.junction_num_nts == args['junction_num_nts']
        ).order_by(desc(Sequence.copy_number)).all()
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
        ).order_by(desc(Sequence.copy_number_in_sample)).all()
        self._print(worker_id, 'Collapsing bucket in subject {} with {} '
                    'seqs'.format(args['subject_id'], len(seqs)))

        collapse_seqs(
            self._session, seqs, 'copy_number_in_sample',
            'copy_number_in_subject', 'collapse_to_subject_seq_id',
            'collapse_to_subject_sample_id'
        )
        if self._tasks % 100 == 0:
            self._session.commit()

    def cleanup(self, worker_id):
        self._print(worker_id, 'Committing collapsed sequences')
        self._session.commit()
        self._session.close()


def seq_to_regex(seq):
    return re.compile(''.join(
        map(lambda c: '[{}N]'.format(c) if c != 'N' else '[ATCGN]', seq)
    ))


def collapse_seqs(session, seqs, copy_field, collapse_copy_field,
                  collapse_seq_id_field, collapse_sample_id_field=None):
    to_process = [{
        'sample_id': s.sample_id,
        'seq_id': s.seq_id,
        'sequence': s.sequence,
        'cn': getattr(s, copy_field)
    } for s in seqs]

    while len(to_process) > 0:
        # Get the largest sequence in the list
        larger = to_process[0]
        # Remove it from the list to process
        to_process = to_process[1:]
        # Compile its sequence into regex
        pattern = seq_to_regex(larger['sequence'])
        # Iterate over all smaller sequences to find matches
        for i in reversed(range(0, len(to_process))):
            smaller = to_process[i]
            if pattern.match(smaller['sequence']) is not None:
                # Add the smaller sequence's copy number to the larger
                larger['cn'] += smaller['cn']
                # If the smaller sequence matches the larger, collapse it to the
                # larger
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

    session.query(Sequence).filter(
        Sequence.sample_id.in_(samples_ids)
    ).update({
        'collapse_to_subject_seq_id': None,
        'collapse_to_subject_sample_id': None,
        'copy_number_in_subject': None,
        'collapse_to_clone_seq_id': None,
        'collapse_to_clone_sample_id': None,
        'copy_number_in_clone': None,
    }, synchronize_session=False)
    session.commit()

    tasks = concurrent.TaskQueue()

    for sample_id in samples_ids:
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
        'copy_number_in_clone': None,
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
                    nproc):
    session = config.init_db(master_db_config, data_db_config)
    for clone_id in clone_ids:
        seqs = session.query(
            Sequence.sample_id,
            Sequence.seq_id,
            Sequence.sequence,
            Sequence.copy_number_in_subject
        ).filter(
            Sequence.clone_id == clone_id,
            Sequence.copy_number_in_subject > 0
        ).order_by(
            desc(Sequence.copy_number_in_subject)
        ).all()

        collapse_seqs(
            session, seqs, 'copy_number_in_subject', 'copy_number_in_clone',
            'collapse_to_clone_seq_id', 'collapse_to_clone_sample_id'
        )

    session.commit()

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
                copy_number_in_subject > 0) as j
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
            s.copy_number_in_clone=0
    '''))
    session.commit()


def run_collapse(session, args):
    mod_log.make_mod('collapse', session=session, commit=True,
                     info=vars(args))
    levels = {
        'sample': (collapse_samples, Sample),
        'subject': (collapse_subjects, Subject),
        'clone': (collapse_clones, Clone)
    }
    ids = args.ids or map(
        lambda r: r.id, session.query(levels[args.level][1]).all()
    )
    print 'Collapsing sequences in {} {}s'.format(
        len(ids), args.level)
    session.close()
    levels[args.level][0](args.master_db_config, args.data_db_config,
                          ids, args.nproc)
