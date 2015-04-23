import sldb.common.config as config
from sldb.common.models import Sample, Sequence
import sldb.util.concurrent as concurrent
from sldb.util.funcs import seq_to_regex

class CollapseWorker(concurrent.Worker):
    def __init__(self, session):
        self._session = session

    def do_task(self, worker_id, sample_id):
        seqs_by_size = {}
        q = self._session.query(
            Sequence.v_gene, Sequence.j_gene, Sequence.junction_num_nts,
            Sequence.seq_id, Sequence.sequence, Sequence.copy_number
        ).filter(
            Sequence.sample_id == sample_id
        ).order_by(Sequence.copy_number)

        for seq in q:
            key = (seq.v_gene, seq.j_gene, seq.junction_num_nts)
            if key not in seqs_by_size:
                seqs_by_size[key] = []
            seqs_by_size[key].append((seq.seq_id, seq.sequence, seq.copy_number))

        self._print(
            worker_id,
            'Collapsing sample {} ({} seqs)'.format(sample_id, sum(map(len,
                seqs_by_size.values()))))
        for i, bucket in enumerate(seqs_by_size.values()):
            self._print(worker_id,
                '   Collapsing sample {} bucket {}/{} ({} seqs)'.format(
                    sample_id, i, len(seqs_by_size), len(bucket)))
            self._collapse_bucket(sample_id, bucket)
        self._print(worker_id, 'Finished entire sample {}'.format(sample_id))
        self._session.commit()

    def _collapse_bucket(self, sample_id, bucket):
        new_cns = {s[0]: s[2] for s in bucket}
        for i, seq1 in enumerate(bucket):
            seq1_seq_id = seq1[0]
            seq1_seq = seq1[1]
            if new_cns[seq1_seq_id] == 0:
                continue
            pattern = seq_to_regex(seq1_seq)
            for j, seq2 in enumerate(bucket[i+1:]):
                seq2_seq_id = seq2[0]
                seq2_seq = seq2[1]
                seq2_copy_number = seq2[2]
                if (new_cns[seq2_seq_id] > 0 
                        and pattern.match(seq2_seq) is not None):
                    new_cns[seq1_seq_id] += seq2_copy_number
                    new_cns[seq2_seq_id] = 0 

        for seq_id, cn in new_cns.iteritems():
            self._session.query(Sequence).filter(
                Sequence.sample_id == sample_id,
                Sequence.seq_id == seq_id
            ).update({
                'copy_number_in_sample': cn
            })


def run_collapse(session, args):
    if args.samples is None:
        sample_ids = map(lambda s: s.id, session.query(Sample.id).all())
    else:
        sample_ids = args.samples

    tasks = concurrent.TaskQueue()
    for sample_id in sample_ids:
        tasks.add_task(sample_id)

    for i in range(0, args.nproc):
        session = config.init_db(args.master_db_config, args.data_db_config)
        tasks.add_worker(CollapseWorker(session))
    tasks.start()
