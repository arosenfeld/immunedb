import argparse
import distance
from collections import Counter

from sqlalchemy import desc, distinct
from sqlalchemy.sql import func

from sldb.identification.identify import VDJSequence
import sldb.util.lookups as lookups
from sldb.util.funcs import page_query
from sldb.common.models import *


def _consensus(strings):
    cons = []
    for chars in zip(*strings):
        cons.append(Counter(chars).most_common(1)[0][0])

    return ''.join(cons)


def _similar_to_all(seq, clone_query, min_similarity):
    for comp_seq in clone_query:
        try:
            if (1 - distance.hamming(comp_seq, seq) / float(len(comp_seq))) < \
                    min_similarity:
                return False
        except:
            return False
    return True


def _assign_identical_sequences(session, replaced_seq, subject_id, clone_id):
    session.query(Sequence)\
        .filter(Sequence.sequence_replaced == replaced_seq,
                Sequence.sample.has(subject_id=subject_id))\
        .update({
            'clone_id': clone_id
        }, synchronize_session=False)


def _get_subject_clones(session, subject_id, min_similarity, limit_alignments,
                        include_indels, per_commit):
    clone_cache = {}
    new_clones = 0
    duplicates = 0
    query = session.query(
        Sequence,
        func.count(Sequence.seq_id).label('others')
        ).filter(
            Sequence.sample.has(subject_id=subject_id),
            Sequence.clone_id.is_(None),
        )
    if not include_indels:
        query = query.filter(
            Sequence.probable_indel_or_misalign == 0)
    query = query.order_by(desc(Sequence.copy_number))

    query = query.group_by(Sequence.sequence_replaced).having(
        func.sum(Sequence.copy_number) > 1)
    for i, seqr in enumerate(query):
        if i > 0 and i % per_commit == 0:
            session.commit()
            print 'Committed {} (new clones={}, duplicates={})'.format(
                i, new_clones, duplicates)

        seq = seqr.Sequence

        if '*' in seq.junction_aa:
            continue
        seq_clone = None
        # Key for cache has implicit subject_id due to function parameter
        key = (seq.v_gene, seq.j_gene, seq.junction_num_nts,
               seq.junction_aa)
        if key in clone_cache:
            seq.clone = clone_cache[key]
            _assign_identical_sequences(session, seq.sequence, subject_id,
                                        clone_cache[key].id)
            continue

        for clone in session.query(Clone)\
                .filter(Clone.subject_id == subject_id,
                        Clone.v_gene == seq.v_gene,
                        Clone.j_gene == seq.j_gene,
                        Clone.cdr3_num_nts == seq.junction_num_nts):
            seqs_in_clone = map(lambda s: s.junction_aa,
                                session.query(Sequence.junction_aa)
                                .filter(Sequence.clone == clone))

            if _similar_to_all(seq.junction_aa, seqs_in_clone,
                               min_similarity):
                seq_clone = clone
                break

        if seq_clone is None:
            new_clone = Clone(subject_id=subject_id,
                              v_gene=seq.v_gene,
                              j_gene=seq.j_gene,
                              cdr3_num_nts=seq.junction_num_nts)
            new_clones += 1
            session.add(new_clone)
            session.flush()
            seq_clone = new_clone

        seq.clone = seq_clone

        if seqr.others > 1:
            duplicates += (seqr.others - 1)
            _assign_identical_sequences(session, seq.sequence_replaced,
                                        subject_id, seq_clone.id)

        clone_cache[key] = seq_clone

    session.commit()


def _assign_clones_to_groups(session, subject_id, per_commit):
    for i, clone in enumerate(session.query(Clone).filter(
            Clone.subject_id == subject_id)):
        seqs = session.query(
            Sequence.junction_nt, Sequence.germline
        ).filter(Sequence.clone_id == clone.id).all()

        clone.cdr3_nt = _consensus(map(lambda s:
                                   s.junction_nt, seqs))
        cdr3_aa = lookups.aas_from_nts(clone.cdr3_nt, '')

        group = session.query(CloneGroup).filter(
            CloneGroup.subject_id == subject_id,
            CloneGroup.v_gene == clone.v_gene,
            CloneGroup.j_gene == clone.j_gene,
            CloneGroup.cdr3_num_nts == clone.cdr3_num_nts,
            CloneGroup.cdr3_aa == cdr3_aa).first()

        if group is None:
            germline = seqs[0].germline

            group = CloneGroup(subject_id=subject_id,
                               v_gene=clone.v_gene,
                               j_gene=clone.j_gene,
                               cdr3_num_nts=clone.cdr3_num_nts,
                               cdr3_aa=cdr3_aa,
                               germline=germline)
            session.add(group)
        clone.group = group

        if i > 0 and i % per_commit == 0:
            print 'Committed {}'.format(i)
            session.commit()

    session.commit()


def run_clones(session, args):
    if args.subjects is None:
        subjects = map(lambda s: s.id, session.query(Subject.id).all())
    else:
        subjects = args.subjects

    for sid in subjects:
        print 'Assigning clones to subject', sid
        _get_subject_clones(session, sid, args.similarity / 100.0,
                            args.limit_alignments, args.include_indels,
                            args.commits)
        print 'Assigning clones to groups'
        _assign_clones_to_groups(session, sid, args.commits)
