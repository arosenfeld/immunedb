import argparse
from collections import Counter
import re

from sqlalchemy import and_, desc, distinct
from sqlalchemy.sql import exists, func, text

import dnautils
from sldb.common.models import (CDR3_OFFSET, Clone, CloneStats, Sequence,
                                SequenceCollapse, Subject)
import sldb.common.modification_log as mod_log
from sldb.identification.identify import VDJSequence
from sldb.identification.v_genes import get_common_seq, VGene, VGermlines
import sldb.util.lookups as lookups


def _consensus(strings):
    """Gets the unweighted consensus from a list of strings

    :param list strings: A set of equal-length strings.

    :returns: A consensus string
    :rtype: str

    """
    chrs = [Counter(chars).most_common(1)[0][0] for chars in zip(*strings)]
    return ''.join(chrs)


def _similar_to_all(seq, rest, min_similarity):
    """Determines if the string ``seq`` is at least ``min_similarity`` similar
    to the list of strings ``rest``.

    :param str seq: The string to compare
    :param list rest: The list of strings to compare to

    :returns: If ``seq`` is similar to every sequence in ``rest``
    :rtype: bool

    """
    for comp_seq in rest:
        sim_frac = 1 - dnautils.hamming(comp_seq, seq) / float(len(comp_seq))
        if sim_frac < min_similarity:
            return False
    return True


def _get_subject_clones(session, subject_id, min_similarity, include_indels,
                        exclude_partials, min_identity, min_copy):
    """Assigns clones to all viable sequences in the specified subject.

    :param Session session: The database session
    :param int subject_id: ID of the subject
    :param float min_similarity: The minimum similarity required for two
        sequences to be in the same clone a clone
    :param bool include_indels: If sequences with possible indels should be
        assigned to clones
    :param bool exclude_partials: If partial sequences should be excluded from
        clonal assignment
    :param float min_identity: The minimum V-gene germline-identity required
        for a sequence to be assigned a clone

    :returns: The list of clone IDs which have been created or updated
    :rtype: list

    """
    clone_cache = {}
    new_clones = 0
    to_update = set([])
    query = session.query(
        Sequence
    ).join(SequenceCollapse).filter(
        Sequence.clone_id.is_(None),

        Sequence.subject_id == subject_id,
        ~Sequence.cdr3_aa.like('%*%'),

        SequenceCollapse.copy_number_in_subject >= min_copy
    )
    if min_identity > 0:
        query = query.filter(
            Sequence.v_match / Sequence.v_length >= min_identity
        )
    if not include_indels:
        query = query.filter(Sequence.probable_indel_or_misalign == 0)
    if exclude_partials:
        query = query.filter(Sequence.partial == 0)

    query = query.order_by(desc(SequenceCollapse.copy_number_in_subject))

    total = query.count()
    for i, seq in enumerate(query):
        if i > 0 and i % 1000 == 0:
            session.commit()
            print 'Committed {}/{} (clones={})'.format(i, total, new_clones)

        # Key for cache has implicit subject_id due to function parameter
        key = (seq.v_gene, seq.j_gene, seq.cdr3_num_nts, seq.cdr3_aa)
        if key in clone_cache:
            seq.clone = clone_cache[key]
            continue

        seq_clone = None
        for clone in session.query(Clone).filter(
                Clone.subject_id == subject_id,
                Clone.v_gene == seq.v_gene,
                Clone.j_gene == seq.j_gene,
                Clone.cdr3_num_nts == seq.cdr3_num_nts):
            seqs_in_clone = session.query(
                Sequence.cdr3_aa
            ).join(SequenceCollapse).filter(
                Sequence.clone == clone,
                SequenceCollapse.copy_number_in_subject >= min_copy
            ).group_by(
                Sequence.cdr3_aa
            )
            seqs_in_clone = map(lambda s: s.cdr3_aa, seqs_in_clone)

            if _similar_to_all(seq.cdr3_aa, seqs_in_clone, min_similarity):
                seq_clone = clone
                break

        if seq_clone is None:
            new_clone = Clone(subject_id=subject_id,
                              v_gene=seq.v_gene,
                              j_gene=seq.j_gene,
                              cdr3_num_nts=seq.cdr3_num_nts)
            new_clones += 1
            session.add(new_clone)
            session.flush()
            seq_clone = new_clone

        seq.clone = seq_clone
        to_update.add(seq_clone.id)

        clone_cache[key] = seq_clone

    session.commit()
    return to_update


def _generate_germline(seqs, clone):
    insertions = set([])
    for seq in seqs:
        if seq.insertions is not None:
            insertions.update(set(seq.insertions))
    clone.insertions = insertions

    for seq in seqs:
        seq.clone_insertions = insertions

    rep_seq = seqs[0]
    rep_ins = rep_seq.insertions or 0
    if rep_ins != 0:
        rep_ins = sum((e[1] for e in rep_ins))
    germline = rep_seq.germline[:CDR3_OFFSET + rep_ins]

    for ins in insertions:
        if ins not in rep_seq.insertions:
            pos, size = ins
            germline = germline[:pos] + ('-' * size) + germline[pos:]
    germline += '-' * clone.cdr3_num_nts

    clone.functional = (
        len(germline) % 3 == 0 and
        not lookups.has_stop(germline)
    )

    j_region = rep_seq.germline.replace('-', '')[-rep_seq.post_cdr3_length:]
    germline += j_region

    return germline


def _generate_consensus(session, subject_id, to_update):
    """Generates consensus CDR3s for clones.

    :param Session session: The database session
    :param int subject_id: The ID of the subject
    :param list to_update: The list of clone IDs to assign to groups

    """
    for i, clone in enumerate(session.query(Clone).filter(
            Clone.id.in_(to_update))):
        seqs = session.query(
            Sequence
        ).join(SequenceCollapse).filter(
            Sequence.clone_id == clone.id,
            SequenceCollapse.copy_number_in_subject > 0
        ).all()

        clone.cdr3_nt = _consensus(map(lambda s: s.cdr3_nt, seqs))
        clone.cdr3_aa = lookups.aas_from_nts(clone.cdr3_nt)
        if i > 0 and i % 1000 == 0:
            print 'Committed {}/{}'.format(i, len(to_update))
            session.commit()

        clone.germline = _generate_germline(seqs, clone)

    session.commit()


def run_clones(session, args):
    """Runs the clone-assignment pipeline stage.

    :param Session session: The database session
    :param Namespace args: The arguments passed to the command

    """
    if args.subject_ids is None:
        subjects = map(lambda s: s.id, session.query(Subject.id).all())
    else:
        subjects = args.subject_ids
    mod_log.make_mod('clones', session=session, commit=True,
                     info=vars(args))

    if args.regen:
        print 'Deleting existing clones'
        session.query(Clone).filter(
            Clone.subject_id.in_(subjects)
        ).delete(synchronize_session=False)
        session.commit()

    for sid in subjects:
        print 'Assigning clones to subject', sid
        to_update = _get_subject_clones(
            session, sid, args.similarity / 100.0,
            args.include_indels, args.exclude_partials,
            args.min_identity / 100.0, args.min_copy)
        if len(to_update) > 0:
            print 'Generating consensus CDR3s'
            _generate_consensus(session, sid, to_update)

    print 'Pushing clone IDs to sample sequences'
    session.connection(mapper=Sequence).execute(text('''
        UPDATE
            sequences AS s
        JOIN sequence_collapse AS c
            ON s.sample_id=c.sample_id AND s.ai=c.seq_ai
        JOIN (SELECT seq_id, ai, clone_id from sequences) as s2
            ON c.collapse_to_subject_seq_ai=s2.ai
        SET s.clone_id=s2.clone_id
        WHERE s.seq_id!=s2.seq_id
    '''))

    session.commit()
