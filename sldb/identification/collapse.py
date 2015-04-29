from sqlalchemy import distinct
from sqlalchemy.sql import desc

import sldb.common.config as config
from sldb.common.models import Sample, Sequence
import sldb.common.modification_log as mod_log
import sldb.util.concurrent as concurrent
import sldb.util.funcs as funcs

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

        funcs.collapse_seqs(
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

        funcs.collapse_seqs(
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


def collapse_samples(master_db_config, data_db_config, samples_to_update,
                     nproc):
    session = config.init_db(master_db_config, data_db_config)
    tasks = concurrent.TaskQueue()
    for sample_id in samples_to_update:
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


def collapse_subjects(master_db_config, data_db_config, samples_to_update,
                      nproc):
    session = config.init_db(master_db_config, data_db_config)
    tasks = concurrent.TaskQueue()

    subject_ids = map(lambda r: r.subject_id, session.query(
        distinct(Sample.subject_id).label('subject_id')
    ).filter(
        Sample.id.in_(samples_to_update)
    ).all())
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


def run_collapse(session, master_db_config, data_db_config, nproc,
                 samples_to_update):
    mod_log.make_mod('collapse', session=session, commit=True,
                     info={
                         'to_update': ','.join(map(str,
                             sorted(samples_to_update))),
                     })

    print 'Collapsing samples...'
    collapse_samples(master_db_config, data_db_config, samples_to_update,
                     nproc)
    print 'Collapsing subjects...'
    collapse_subjects(master_db_config, data_db_config, samples_to_update,
                      nproc)
