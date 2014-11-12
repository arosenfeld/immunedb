import argparse
from collections import Counter

import sldb.util.lookups as lookups
from sldb.common.models import *


def _page_query(q, per=10000):
    offset = 0
    while True:
        r = False
        for elem in q.limit(per).offset(offset):
           r = True
           yield elem
        offset += per
        if not r:
            break


def _consensus(strings):
    cons = []
    for chars in zip(*strings):
        cons.append(Counter(chars).most_common(1)[0][0])

    return ''.join(cons)


def _similar_to_all(seq, cluster_query, min_similarity):
    for comp_seq in cluster_query:
        if distance.hamming(cid_seq.sequence, seq.sequence) < min_similarity:
            return False
    return True

def _cluster_subject(session, subject_id, min_similarity):
    print Sample.__table_args__
    for seq in _page_query(session.query(Sequence.sequence)\
            .filter(Sequence.sample.has(subject_id=subject_id),
                    Sequence.cluster_id.is_(None))):
        # NOTE: Would a join be faster here?
        for cluster in session.query(Cluster)\
                .filter(Cluster.subject_id == subject_id,
                        Cluster.v_gene == seq.v_call,
                        Cluster.j_gene == seq.j_call,
                        Cluster.cdr3_num_nts == len(seq.junction_nt)):
            seqs_in_cluster = session.query(Sequence.sequence).filter(
                                   Sequence.cluster == cluster)
            if _similar_to_all(seq, seqs_in_cluster):
                seq.cluster = cluster
                break

        if seq.cluster is None:
            new_cluster = Cluster(subject_id=subject_id,
                                  v_gene=seq.v_call,
                                  j_gene=seq.j_call,
                                  cdr3_num_nts=len(seq.junction_nt))
            session.add(new_cluster)
            seq.cluster = new_cluster


def _assign_clones_to_subject(session, subject_id, per_commit):
    for i, cluster in enumerate(session.query(Cluster).filter(
            Cluster.subject.has(id=subject_id))):
        seqs = session.query(Sequence).filter(
            Sequence.sample.has(subject_id=subject_id),
            Sequence.cluster == cluster).all()

        cdr3_nt = _consensus(map(lambda e: e.junction_nt, seqs))
        cdr3_aa = lookups.aas_from_nts(clone.cdr3_nt, '')
        cluster.cdr3_nt = cdr3_nt

        clone = session.query(Clone).filter(
            Clone.subject_id == subject_id,
            Clone.v_gene == Cluster.v_gene,
            Clone.j_gene == Cluster.j_gene,
            Clone.cdr3_num_nts == Cluster.cdr3_num_nts,
            Clone.cdr3_aa == cdr3_aa).first()

        if clone is None:
            clone = Clone(v_gene=Cluster.v_gene,
                          j_gene=Cluster.j_gene,
                          cdr3_num_nts=Cluster.cdr3_num_nts,
                          cdr3_aa=cdr3_aa)
            session.add(clone)
        cluster.clone = clone

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
        _cluster_subject(session, sid, args.similarity)
        _assign_clones_to_subject(session, sid, args.commits)
