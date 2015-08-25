from sqlalchemy import func

import sldb.util.funcs as funcs
from sldb.identification import AlignmentException, SequenceRecord
from sldb.identification.j_genes import JGermlines
from sldb.identification.v_genes import VGermlines
from sldb.identification.vdj_sequence import VDJSequence
from sldb.common.models import (HashExtension, DuplicateSequence, NoResult,
                                Sequence)


def run_fix_sequences(session, args):
    v_germlines = VGermlines(args.v_germlines)
    j_germlines = JGermlines(args.j_germlines, args.upstream_of_cdr3,
                             args.anchor_len, args.min_anchor_len)

    mutation_cache = {}
    indels = session.query(Sequence).filter(
        Sequence.probable_indel_or_misalign == 1)
    total = indels.count()
    fixed = 0
    print 'Locally aligning {} indels'.format(total)
    for i, seq in enumerate(funcs.periodic_commit(session, indels)):
        try:
            if seq.sample_id not in mutation_cache:
                mutation_cache[seq.sample_id] = session.query(
                    func.avg(Sequence.v_mutation_fraction),
                    func.avg(Sequence.v_length)
                ).filter(Sequence.sample == seq.sample).first()
            avg_mut, avg_len = mutation_cache[seq.sample_id]

            v = VDJSequence(
                seq.seq_id, seq.original_sequence, v_germlines,
                j_germlines, quality=funcs.quality_to_ord(
                    seq.original_quality
                ), locally_align=(avg_mut, avg_len))
            if not v.has_possible_indel:
                existing = session.query(
                    Sequence.seq_id,
                    Sequence.sample_id
                ).filter(
                    Sequence.sample_seq_hash == HashExtension.hash_fields(
                        (seq.sample_id, v.sequence))
                ).first()
                if existing is not None:
                    existing.copy_number += seq.copy_number
                    session.update(DuplicateSequence).filter(
                        DuplicateSequence.sample_id == seq.sample_id,
                        DuplicateSequence.duplicate_seq_id == seq.seq_id
                    ).update({
                        'duplicate_seq_id': existing.seq_id
                    })
                    session.add(DuplicateSequence(
                        duplicate_seq_id=existing.seq_id,
                        sample_id=existing.sample_id,
                        seq_id=seq.seq_id))
                    session.delete(seq)
                else:
                    r = SequenceRecord(v.sequence, seq.sample)
                    r.seq_ids = map(lambda e: e.seq_id, session.query(
                        DuplicateSequence.seq_id
                    ).filter(
                        DuplicateSequence.sample_id == seq.sample_id,
                        DuplicateSequence.duplicate_seq_id == seq.seq_id
                    ).all())
                    r.seq_ids.append(v.id)
                    r.vdj = v
                    session.delete(seq)
                    session.flush()
                    r.add_as_sequence(session, seq.sample, seq.paired)
                fixed += 1
        except AlignmentException as e:
            pass
        if i > 0 and i % 10 == 0:
            print ('Fixed {} of {} processed indel sequences, {} '
                   'remaining').format(fixed, i + 1, total - i + 1)

    session.commit()
