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
from sldb.util.funcs import page_query, seq_to_regex
import sldb.util.lookups as lookups


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
                        include_indels, min_identity, order):
    clone_cache = {}
    new_clones = 0
    to_update = set([])
    duplicates = 0
    query = session.query(
        Sequence,
        func.count(Sequence.seq_id).label('others')
        ).filter(
            Sequence.sample.has(subject_id=subject_id),
            Sequence.v_match / Sequence.v_length >= min_identity,
            Sequence.alignment.in_(limit_alignments),
            Sequence.clone_id.is_(None),
        )
    if not include_indels:
        query = query.filter(
            Sequence.probable_indel_or_misalign == 0)
    if order:
        query = query.order_by(desc(Sequence.copy_number))

    query = query.group_by(Sequence.sequence_replaced).having(
        func.sum(Sequence.copy_number) > 1)

    for i, seqr in enumerate(query):
        if i > 0 and i % 1000 == 0:
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
            _assign_identical_sequences(session, seq.sequence_replaced,
                                        subject_id, clone_cache[key].id)
            continue

        for clone in session.query(Clone)\
                .filter(Clone.subject_id == subject_id,
                        Clone.v_gene == seq.v_gene,
                        Clone.j_gene == seq.j_gene,
                        Clone.cdr3_num_nts == seq.junction_num_nts):
            seqs_in_clone = map(lambda s: s.junction_aa,
                                session.query(Sequence.junction_aa)
                                .filter(Sequence.clone == clone).group_by(
                                    Sequence.junction_aa
                                ))

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
        to_update.add(seq_clone.id)

        if seqr.others > 1:
            duplicates += (seqr.others - 1)
            _assign_identical_sequences(session, seq.sequence_replaced,
                                        subject_id, seq_clone.id)

        clone_cache[key] = seq_clone

    session.commit()
    return to_update


def _assign_clones_to_groups(session, subject_id, to_update):
    if len(to_update) == 0:
        return
    for i, clone in enumerate(session.query(Clone).filter(
            Clone.id.in_(to_update))):
        seqs = session.query(Sequence).filter(
            Sequence.clone_id == clone.id).all()

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
        session.query(Sequence).filter(
            Sequence.clone_id == clone_id
        ).update({
            'copy_number_in_clone': 0
        })
        session.flush()

        seqs_by_size = session.query(
            Sequence.seq_id, Sequence.sample_id, Sequence.sequence,
            Sequence.copy_number_in_sample
        ).filter(
            Sequence.clone_id == clone_id,
            Sequence.copy_number_in_sample > 0
        ).order_by(Sequence.copy_number_in_sample).all()
        new_cns = {(s.sample_id, s.seq_id): s.copy_number_in_sample for s in
            seqs_by_size}

        for i, seq1 in enumerate(seqs_by_size):
            s1_key = (seq1.sample_id, seq1.seq_id)
            if new_cns[s1_key] == 0:
                continue
            pattern = seq_to_regex(seq1.sequence)
            for j, seq2 in enumerate(seqs_by_size[i+1:]):
                s2_key = (seq2.sample_id, seq2.seq_id)
                if (new_cns[s2_key] > 0 
                        and pattern.match(seq2.sequence) is not None):
                    new_cns[s1_key] += seq2.copy_number_in_sample
                    new_cns[s2_key] = 0 

        for (sample_id, seq_id), cn in new_cns.iteritems():
            session.query(Sequence).filter(
                Sequence.sample_id == sample_id,
                Sequence.seq_id == seq_id
            ).update({
                'copy_number_in_clone': cn
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
            args.min_identity / 100.0, args.order)
        print 'Assigning clones to groups'
        _assign_clones_to_groups(session, sid, to_update)
        _collapse_sequences(session, to_update)
