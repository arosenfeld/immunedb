import argparse
from collections import Counter
import re

import distance

from sqlalchemy import desc, distinct
from sqlalchemy.sql import func

from sldb.common.models import *
import sldb.common.modification_log as mod_log
from sldb.identification.identify import VDJSequence
from sldb.identification.v_genes import VGene
from sldb.util.funcs import collapse_seqs, page_query, seq_to_regex
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
        Sequence.v_match / Sequence.v_length >= min_identity,
        Sequence.alignment.in_(limit_alignments),
        ~Sequence.junction_aa.like('%*%'),

        Sequence.copy_number_in_subject > 1
    )
    if not include_indels:
        query = query.filter(Sequence.probable_indel_or_misalign == 0)

    query = query.order_by(desc(Sequence.copy_number_in_subject))

    for i, seq in enumerate(query):
        if i > 0 and i % 1000 == 0:
            session.commit()
            print 'Committed {} (new clones={})'.format( i, new_clones)

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
                Sequence.copy_number_in_subject > 0
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
            print 'Committed {}'.format(i)
            session.commit()

    session.commit()


def _collapse_sequences(session, to_update):
    print 'Collapsing sequences'
    for clone_id in to_update:
        seqs = session.query(
            Sequence.sample_id,
            Sequence.seq_id,
            Sequence.sequence,
            Sequence.copy_number_in_subject
        ).filter(
            Sequence.clone_id == clone_id,
        ).order_by(
            desc(Sequence.copy_number_in_subject)
        ).all()

        collapse_seqs(
            session, seqs, 'copy_number_in_subject', 'copy_number_in_clone',
            'collapse_to_clone_seq_id', 'collapse_to_clone_sample_id'
        )

    session.commit()

def _push_clones_down(session, to_update):
    clone_assigned = session.query(
        Sequence.sample_id,
        Sequence.seq_id,
        Sequence.clone_id,
        Sequence.copy_number,
    ).filter(
        Sequence.copy_number_in_clone > 0,
        Sequence.clone_id.in_(to_update)
    )

    for clone_seq in clone_assigned:
        subject_assigned = session.query(
            Sequence
        ).filter(
            Sequence.collapse_to_clone_sample_id == clone_seq.sample_id,
            Sequence.collapse_to_clone_seq_id == clone_seq.seq_id
        )
        for subject_seq in subject_assigned:
            subject_seq.clone_id = clone_seq.clone_id
            sample_assigned = session.query(
                Sequence
            ).filter(
                Sequence.collapse_to_subject_sample_id == subject_seq.sample_id,
                Sequence.collapse_to_subject_seq_id == subject_seq.seq_id
            )
            for sample_seq in sample_assigned:
                # Update the sequences and everything collapsed to it in its sample
                sample_seq.clone_id = clone_seq.clone_id
                session.query(Sequence).filter(
                    Sequence.collapse_to_sample_seq_id == sample_seq.seq_id,
                    Sequence.sample_id == sample_seq.sample_id
                ).update({
                    'clone_id': clone_seq.clone_id
                })
    session.commit()

def run_clones(session, args):
    if args.subjects is None:
        subjects = map(lambda s: s.id, session.query(Subject.id).all())
    else:
        subjects = args.subjects
    mod_log.make_mod('clones', session=session, commit=True,
                     info=vars(args))

    for sid in subjects:
        print 'Assigning clones to subject', sid
        to_update = _get_subject_clones(
            session, sid, args.similarity / 100.0,
            args.limit_alignments, args.include_indels,
            args.min_identity / 100.0)
        print 'Assigning clones to groups'
        _assign_clones_to_groups(session, sid, to_update)
        to_update = map(lambda r:r.id, session.query(Clone.id).filter(
            Clone.id != None).all())
        _collapse_sequences(session, to_update)
        print 'Pushing clone IDs down'
        _push_clones_down(session, to_update)
