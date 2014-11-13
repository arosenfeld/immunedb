import argparse
import distance
from collections import Counter

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
        if distance.hamming(comp_seq, seq) / float(len(comp_seq)) \
                < min_similarity:
            return False
    return True

def _get_subject_clones(session, subject_id, min_similarity, per_commit):
    clone_cache = {}
    for i, seq in enumerate(session.query(Sequence)\
            .filter(Sequence.sample.has(subject_id=subject_id),
                    Sequence.copy_number_iden > 1,
                    Sequence.clone_id.is_(None))):

        # Key for cache has implicit subject_id due to function parameter
        key = (seq.v_call, seq.j_call, seq.junction_aa)
        if key in clone_cache:
            seq.clone = clone_cache[key]
            continue

        for clone in session.query(Clone)\
                .filter(Clone.subject_id == subject_id,
                        Clone.v_gene == seq.v_call,
                        Clone.j_gene == seq.j_call,
                        Clone.cdr3_num_nts == len(seq.junction_nt)):
            seqs_in_clone = map(lambda s: s.sequence,
                                  session.query(Sequence.sequence).filter(
                                    Sequence.clone == clone))
            if _similar_to_all(seq.sequence, seqs_in_clone, min_similarity):
                seq.clone = clone
                break

        if seq.clone is None:
            new_clone = Clone(subject_id=subject_id,
                                  v_gene=seq.v_call,
                                  j_gene=seq.j_call,
                                  cdr3_num_nts=len(seq.junction_nt))
            session.add(new_clone)
            seq.clone = new_clone

        clone_cache[key] = seq.clone

        if i > 0 and i % per_commit == 0:
            session.commit()
            print 'Committed {}'.format(i)

    session.commit()

def _assign_clones_to_groups(session, subject_id, per_commit):
    for i, clone in enumerate(session.query(Clone).filter(
            Clone.subject_id == subject_id)):
        seqs = session.query(Sequence).filter(
            Sequence.clone == clone).all()

        clone.cdr3_nt = _consensus(map(lambda e: e.junction_nt, seqs))
        cdr3_aa = lookups.aas_from_nts(clone.cdr3_nt, '')

        group = session.query(CloneGroup).filter(
            CloneGroup.subject_id == subject_id,
            CloneGroup.v_gene == clone.v_gene,
            CloneGroup.j_gene == clone.j_gene,
            CloneGroup.cdr3_num_nts == clone.cdr3_num_nts,
            CloneGroup.cdr3_aa == cdr3_aa).first()

        if group is None:
            group = CloneGroup(subject_id=subject_id,
                               v_gene=clone.v_gene,
                               j_gene=clone.j_gene,
                               cdr3_num_nts=clone.cdr3_num_nts,
                               cdr3_aa=cdr3_aa)
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
        _get_subject_clones(session, sid, args.similarity / 100.0, args.commits)
        print 'Assigning clones to groups'
        _assign_clones_to_groups(session, sid, args.commits)
