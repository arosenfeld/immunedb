import argparse
import distance
from collections import Counter

from sqlalchemy import distinct
from sqlalchemy.sql import func

import sldb.identification.germlines as germlines
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
        if distance.hamming(comp_seq, seq) / float(len(comp_seq)) \
                < min_similarity:
            return False
    return True

def _get_subject_clones(session, subject_id, min_similarity, per_commit):
    clone_cache = {}
    new_clones = 0
    samples = map(lambda s: s.id, session.query(Sample).filter(
        Sample.subject_id == subject_id))
    for i, seq in enumerate(session.query(SequenceMapping).filter(
                SequenceMapping.sample_id.in_(samples),
                SequenceMapping.clone_id.is_(None),
                SequenceMapping.copy_number > 1)):
        if i > 0 and i % per_commit == 0:
            session.commit()
            print 'Committed {} (new clones={})'.format(i, new_clones)
        seq_iden = seq.identity_seq
        if '*' in seq_iden.junction_aa:
            continue
        seq_clone = None
        # Key for cache has implicit subject_id due to function parameter
        key = (seq_iden.v_call, seq_iden.j_call, len(seq_iden.junction_nt), 
               seq_iden.junction_aa)
        if key in clone_cache:
            seq.clone = clone_cache[key]
            continue

        for clone in session.query(Clone)\
                .filter(Clone.subject_id == subject_id,
                        Clone.v_gene == seq_iden.v_call,
                        Clone.j_gene == seq_iden.j_call,
                        Clone.cdr3_num_nts == len(seq_iden.junction_nt)):
            seqs_in_clone = map(lambda s: s.identity_seq.junction_aa,
                                session.query(SequenceMapping)\
                                    .filter(SequenceMapping.clone == clone))

            if _similar_to_all(seq_iden.junction_aa, seqs_in_clone, min_similarity):
                seq_clone = clone
                break

        if seq_clone is None:
            new_clone = Clone(subject_id=subject_id,
                                  v_gene=seq_iden.v_call,
                                  j_gene=seq_iden.j_call,
                                  cdr3_num_nts=len(seq_iden.junction_nt))
            new_clones += 1
            session.add(new_clone)
            session.flush()
            seq_clone = new_clone

        session.query(SequenceMapping)\
            .filter(SequenceMapping.identity_seq_id == seq.identity_seq_id,
                    SequenceMapping.sample_id.in_(samples),
                    SequenceMapping.copy_number > 1)\
            .update({
                'clone_id': seq_clone.id
            }, synchronize_session='fetch')

        clone_cache[key] = seq_clone

    session.commit()

def _assign_clones_to_groups(session, subject_id, per_commit):
    for i, clone in enumerate(session.query(Clone).filter(
            Clone.subject_id == subject_id)):
        seqs = session.query(SequenceMapping).filter(
            SequenceMapping.clone_id == clone.id).all()

        clone.cdr3_nt = _consensus(map(lambda e: 
                                   e.identity_seq.junction_nt, seqs))
        cdr3_aa = lookups.aas_from_nts(clone.cdr3_nt, '')

        group = session.query(CloneGroup).filter(
            CloneGroup.subject_id == subject_id,
            CloneGroup.v_gene == clone.v_gene,
            CloneGroup.j_gene == clone.j_gene,
            CloneGroup.cdr3_num_nts == clone.cdr3_num_nts,
            CloneGroup.cdr3_aa == cdr3_aa).first()

        if group is None:
            # NOTE: This is for v-ties, and there may be a better way
            v = clone.v_gene.split('|')[0]
            germline = '{}{}{}'.format(
                germlines.v[v][0:VDJSequence.CDR3_OFFSET],
                '-' * clone.cdr3_num_nts,
                germlines.j[clone.j_gene])
            germline = germline[:len(seqs[0].identity_seq.sequence_replaced)]

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
        #print 'Assigning clones to subject', sid
        #_get_subject_clones(session, sid, args.similarity / 100.0, args.commits)
        print 'Assigning clones to groups'
        _assign_clones_to_groups(session, sid, args.commits)
