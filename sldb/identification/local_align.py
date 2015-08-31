import multiprocessing as mp
import Queue
from sqlalchemy import func

import sldb.util.funcs as funcs
from sldb.identification import AlignmentException, SequenceRecord
from sldb.identification.j_genes import JGermlines
from sldb.identification.v_genes import VGermlines
from sldb.identification.vdj_sequence import VDJSequence
from sldb.common.models import (HashExtension, DuplicateSequence, NoResult,
                                Sample, Sequence)
import sldb.util.concurrent as concurrent


class LocalAlignmentWorker(concurrent.Worker):
    def __init__(self, v_germlines, j_germlines, completed_tasks):
        self._v_germlines = v_germlines
        self._j_germlines = j_germlines
        self._completed_tasks = completed_tasks

    def do_task(self, args):
        seq = args['seq']
        try:
            v = VDJSequence(
                seq.seq_id, seq.original_sequence, self._v_germlines,
                self._j_germlines, quality=funcs.quality_to_ord(
                    seq.original_quality
                ), locally_align=(args['avg_mut'], args['avg_len'])
            )

            if not v.has_possible_indel:
                self._print('Adding {}'.format(v.id))
                args.update({
                    'vdj': v
                })
                self._completed_tasks.put(args)
        except AlignmentException:
            pass


def run_fix_sequences(session, args):
    v_germlines = VGermlines(args.v_germlines)
    j_germlines = JGermlines(args.j_germlines, args.upstream_of_cdr3,
                             args.anchor_len, args.min_anchor_len)

    mutation_cache = {}
    tasks = concurrent.TaskQueue()
    completed_tasks = mp.Queue()

    indels = session.query(Sequence).filter(
        Sequence.probable_indel_or_misalign == 1)
    total = indels.count()
    fixed = 0
    print 'Creating task queue for {} indels'.format(total)
    for i, seq in enumerate(indels):
        if seq.sample_id not in mutation_cache:
            mutation_cache[seq.sample_id] = session.query(
                func.avg(Sequence.v_mutation_fraction),
                func.avg(Sequence.v_length)
            ).filter(Sequence.sample == seq.sample).first()
        avg_mut, avg_len = mutation_cache[seq.sample_id]
        session.expunge(seq)
        tasks.add_task({
            'sample_id': seq.sample_id,
            'seq': seq,
            'copy_number': seq.copy_number,
            'paired': seq.paired,
            'avg_mut': avg_mut,
            'avg_len': avg_len
        })

    for i in range(0, min(args.nproc, tasks.num_tasks)):
        tasks.add_worker(LocalAlignmentWorker(v_germlines, j_germlines,
                                              completed_tasks))

    tasks.start()
    i = 0
    while True:
        i += 1
        try:
            task = completed_tasks.get()
        except Queue.Empty:
            break

        vdj = task['vdj']
        sample_id = task['sample_id']
        copy_number = task['copy_number']
        paired = task['paired']

        # Delete the sequence which was successfully locally aligned
        session.query(Sequence).filter(
            Sequence.sample_id == sample_id,
            Sequence.seq_id == vdj.id
        ).delete()
        session.flush()

        # Check if there is an existing sequence with the same aligned sequence
        existing = session.query(
            Sequence.seq_id,
            Sequence.sample_id
        ).filter(
            Sequence.sample_seq_hash == HashExtension.hash_fields(
                (sample_id, vdj.sequence))
        ).first()
        if existing is not None:
            # There is an existing sequence; add the aligned sequence as a
            # duplicate of it.
            existing.copy_number += copy_number
            # Set duplicates of the aligned sequence to duplicates of the
            # existing sequence
            session.update(DuplicateSequence).filter(
                DuplicateSequence.sample_id == sample_id,
                DuplicateSequence.duplicate_seq_id == seq_id
            ).update({
                'duplicate_seq_id': existing.seq_id
            })

            # Add the sequence as a duplicate
            session.add(DuplicateSequence(
                duplicate_seq_id=existing.seq_id,
                sample_id=existing.sample_id,
                seq_id=vdj.id))
        else:
            # The aligned sequence is unique; add it as a regular sequence
            r = SequenceRecord(vdj.sequence, sample_id)
            r.seq_ids = map(lambda e: e.seq_id, session.query(
                DuplicateSequence.seq_id
            ).filter(
                DuplicateSequence.sample_id == sample_id,
                DuplicateSequence.duplicate_seq_id == vdj.id
            ).all())
            r.seq_ids.append(vdj.id)
            r.vdj = vdj
            sample = session.query(Sample).filter(Sample.id == sample_id).one()
            r.add_as_sequence(session, sample, paired)
        fixed += 1

        if i > 0 and i % 10 == 0:
            print 'Processed {} indels'.format(i)
            session.commit()
