import argparse
from collections import Counter
import re

import distance

from sqlalchemy import and_, desc, distinct
from sqlalchemy.sql import exists, func

from sldb.common.models import *
import sldb.common.modification_log as mod_log
from sldb.identification.identify import VDJSequence
from sldb.identification.v_genes import VGene
import sldb.util.lookups as lookups


def _consensus(strings):
    chrs = [Counter(chars).most_common(1)[0][0] for chars in zip(*strings)]
    return ''.join(chrs)


def _similar_to_all(seq, clone_query, min_similarity):
    for comp_seq in clone_query:
        sim_frac = 1 - distance.hamming(comp_seq, seq) / float(len(comp_seq))
        if sim_frac < min_similarity:
            return False
    return True


def _get_subject_clones(session, subject_id, min_similarity, limit_alignments,
                        include_indels, min_identity):
    clone_cache = {}
    new_clones = 0
    to_update = set([])
    query = session.query(
        Sequence
    ).filter(
        Sequence.clone_id.is_(None),

        Sequence.sample.has(subject_id=subject_id),
        Sequence.alignment.in_(limit_alignments),
        ~Sequence.junction_aa.like('%*%'),

        Sequence.copy_number_in_subject > 1
    )
    if min_identity > 0:
        query = query.filter(
            Sequence.v_match / Sequence.v_length >= min_identity
        )
    if not include_indels:
        query = query.filter(Sequence.probable_indel_or_misalign == 0)

    query = query.order_by(desc(Sequence.copy_number_in_subject))

    total = query.count()
    for i, seq in enumerate(query):
        if i > 0 and i % 1000 == 0:
            session.commit()
            print 'Committed {}/{} (new clones={})'.format(i, total,
                                                           new_clones)

        # Key for cache has implicit subject_id due to function parameter
        key = (seq.v_gene, seq.j_gene, seq.junction_num_nts,
               seq.junction_aa)
        if key in clone_cache:
            seq.clone = clone_cache[key]
            continue

        seq_clone = None
        for clone in session.query(Clone)\
                .filter(Clone.subject_id == subject_id,
                        Clone.v_gene == seq.v_gene,
                        Clone.j_gene == seq.j_gene,
                        Clone.cdr3_num_nts == seq.junction_num_nts):
            seqs_in_clone = session.query(
                Sequence.junction_aa
            ).filter(
                Sequence.clone == clone,
                Sequence.copy_number_in_subject > 1
            ).group_by(
                Sequence.junction_aa
            )
            seqs_in_clone = map(lambda s: s.junction_aa, seqs_in_clone)

            if _similar_to_all(seq.junction_aa, seqs_in_clone, min_similarity):
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
        to_update.add(seq_clone.id)

        clone_cache[key] = seq_clone

    session.commit()
    return to_update


def _assign_clones_to_groups(session, subject_id, to_update):
    if len(to_update) == 0:
        return
    for i, clone in enumerate(session.query(Clone).filter(
            Clone.id.in_(to_update))):
        seqs = session.query(Sequence).filter(
            Sequence.clone_id == clone.id,
            Sequence.copy_number_in_subject > 0
        ).all()

        clone.cdr3_nt = _consensus(map(lambda s:
                                   s.junction_nt, seqs))
        cdr3_aa = lookups.aas_from_nts(clone.cdr3_nt)

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

        if i > 0 and i % 1000 == 0:
            print 'Committed {}/{}'.format(i, len(to_update))
            session.commit()

    session.commit()


def run_clones(session, args):
    if args.subject_ids is None:
        subjects = map(lambda s: s.id, session.query(Subject.id).all())
    else:
        subjects = args.subject_ids
    mod_log.make_mod('clones', session=session, commit=True,
                     info=vars(args))

    if args.regen:
        print 'Unassigning existing clones'
        session.query(Sequence).filter(
            exists().where(
                and_(Clone.id == Sequence.clone_id,
                     Clone.subject_id.in_(subjects))
            )
        ).update({
            'clone_id': None,
            'mutations_from_clone': None,
            'collapse_to_clone_seq_id': None,
            'collapse_to_clone_sample_id': None,
            'copy_number_in_clone': None
        }, synchronize_session=False)
        print 'Deleting existing clone stats'
        session.query(CloneStats).filter(
            exists().where(
                and_(CloneStats.clone_id == Clone.id,
                     Clone.subject_id.in_(subjects))
            )
        ).delete(synchronize_session=False)
        print 'Deleting existing clones'
        session.query(Clone).filter(
            Clone.subject_id.in_(subjects)
        ).delete(synchronize_session=False)
        session.commit()

    for sid in subjects:
        print 'Assigning clones to subject', sid
        to_update = _get_subject_clones(
            session, sid, args.similarity / 100.0,
            args.limit_alignments, args.include_indels,
            args.min_identity / 100.0)
        print 'Assigning clones to groups'
        _assign_clones_to_groups(session, sid, to_update)
