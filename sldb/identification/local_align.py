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
import sldb.util.lookups as lookups
import sldb.util.concurrent as concurrent


class LocalAlignmentWorker(concurrent.Worker):
    def __init__(self, v_germlines, j_germlines, completed_tasks,
                 max_deletions, max_insertions):
        self._v_germlines = v_germlines
        self._j_germlines = j_germlines
        self._max_deletions = max_deletions
        self._max_insertions = max_insertions
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

            if (not v.has_possible_indel and len(v.insertions) <=
                    self._max_insertions and len(v.deletions) <=
                    self._max_deletions):
                args.update({
                    'vdj': v
                })
                self._completed_tasks.put(args)
        except AlignmentException:
            pass

    def cleanup(self):
        self._completed_tasks.put(None)


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

    workers = min(args.nproc, tasks.num_tasks)
    for i in range(0, workers):
        tasks.add_worker(LocalAlignmentWorker(v_germlines, j_germlines,
                                              completed_tasks,
                                              args.max_deletions,
                                              args.max_insertions))

    tasks.start(block=False)
    i = 0
    stops = 0
    while True:
        i += 1
        try:
            task = completed_tasks.get()
        except Queue.Empty:
            break
        if task is None:
            stops += 1
            if stops >= workers:
                break
            continue

        vdj = task['vdj']
        sample_id = task['sample_id']
        copy_number = task['copy_number']
        paired = task['paired']

        # Check if there is an existing sequence with the same aligned sequence
        existing = session.query(Sequence).filter(
            Sequence.sample_seq_hash == HashExtension.hash_fields(
                (sample_id, vdj.sequence)),
            Sequence.seq_id != vdj.id
        ).first()
        if existing is not None:
            # There is an existing sequence; add the aligned sequence as a
            # duplicate of it.
            existing.copy_number += copy_number
            # Set duplicates of the aligned sequence to duplicates of the
            # existing sequence
            for dup in session.query(DuplicateSequence).filter(
                    DuplicateSequence.sample_id == sample_id,
                    DuplicateSequence.duplicate_seq_id == vdj.id):
                        dup.duplicate_seq_id, existing.seq_id)
                dup.duplicate_seq_id = existing.seq_id

            # Add the sequence as a duplicate
            session.add(DuplicateSequence(
                duplicate_seq_id=existing.seq_id,
                sample_id=sample_id,
                seq_id=vdj.id))

            # Delete the original sequence
            session.query(Sequence).filter(
                Sequence.sample_id == sample_id,
                Sequence.seq_id == vdj.id
            ).delete()
        else:
            # The aligned sequence is unique; add it as a regular sequence
            old_seq = session.query(Sequence).filter(
                Sequence.seq_id == vdj.id,
                Sequence.sample_id == sample_id
            ).one()

            old_seq.probable_indel_or_misalign = vdj.has_possible_indel
            old_seq.deletions = vdj.deletions
            old_seq.insertions = vdj.insertions
            old_seq.v_gene = funcs.format_ties(vdj.v_gene, 'IGHV')
            old_seq.j_gene = funcs.format_ties(vdj.j_gene, 'IGHJ')

            old_seq.num_gaps = vdj.num_gaps
            old_seq.pad_length = vdj.pad_length

            old_seq.v_match = vdj.v_match
            old_seq.v_length = vdj.v_length
            old_seq.j_match = vdj.j_match
            old_seq.j_length = vdj.j_length

            old_seq.removed_prefix = vdj.removed_prefix
            old_seq.removed_prefix_qual = funcs.ord_to_quality(
               vdj.removed_prefix_qual)
            old_seq.v_mutation_fraction = vdj.mutation_fraction

            old_seq.pre_cdr3_length = vdj.pre_cdr3_length
            old_seq.pre_cdr3_match = vdj.pre_cdr3_match
            old_seq.post_cdr3_length = vdj.post_cdr3_length
            old_seq.post_cdr3_match = vdj.post_cdr3_match

            old_seq.in_frame = vdj.in_frame
            old_seq.functional = vdj.functional
            old_seq.stop = vdj.stop

            old_seq.cdr3_nt = vdj.cdr3
            old_seq.cdr3_num_nts = len(vdj.cdr3)
            old_seq.cdr3_aa = lookups.aas_from_nts(vdj.cdr3)

            old_seq.sequence = str(vdj.sequence)
            old_seq.quality = funcs.ord_to_quality(vdj.quality)

            old_seq.germline = vdj.germline

        fixed += 1

        if i > 0 and i % 10 == 0:
            print 'Processed {} indels'.format(i)
            session.commit()
    session.commit()
