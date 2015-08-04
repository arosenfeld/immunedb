import hashlib
import re

from sldb.common.models import Sequence


def periodic_commit(session, query, interval=100):
    for i, r in enumerate(query):
        if i > 0 and i % interval == 0:
            session.commit()
        yield r
    session.commit()


def trace_seq_collapses(session, seq):
    ret = {}
    sample_col = session.query(
        Sequence.seq_id,
        Sequence.sample_id,
        Sequence.copy_number_in_sample,

        Sequence.collapse_to_subject_seq_id,
        Sequence.collapse_to_subject_sample_id,
    ).filter(
        Sequence.seq_id == seq.collapse_to_sample_seq_id,
        Sequence.sample_id == seq.sample_id
    ).first()

    if sample_col is not None:
        ret['sample'] = {
            'seq_id': sample_col.seq_id,
            'sample_id': sample_col.sample_id,
            'copy_number': sample_col.copy_number_in_sample
        }

        subject_col = session.query(
            Sequence.seq_id,
            Sequence.sample_id,
            Sequence.copy_number_in_subject
        ).filter(
            Sequence.seq_id == sample_col.collapse_to_subject_seq_id,
            Sequence.sample_id == sample_col.collapse_to_subject_sample_id,
        ).first()

        if subject_col is not None:
            ret['subject'] = {
                'seq_id': subject_col.seq_id,
                'sample_id': subject_col.sample_id,
                'copy_number': subject_col.copy_number_in_subject
            }

    return ret


def get_or_create(session, model, **kwargs):
    """Gets or creates a record based on some kwargs search parameters"""
    instance = session.query(model).filter_by(**kwargs).first()
    if instance:
        return instance, False
    else:
        instance = model(**kwargs)
        session.add(instance)
        return instance, True


def find_streak_position(s1, s2, max_streak):
    '''Finds the first streak of max_streak characters where s1 does not equal
    s2

    For example if max_streak is 3:
        ATCGATCGATCGATCG
        ATCGATCGATCTTACG
                     ^--- Returned index
    '''
    streak = 0
    for i, (c1, c2) in enumerate(zip(s1, s2)):
        streak = streak + 1 if c1 != c2 else 0
        if streak >= max_streak:
            return i
    return None


def format_ties(ties, name):
    if ties is None:
        return None
    ties = map(lambda e: e.replace(name, ''), ties)
    return '{}{}'.format(name, '|'.join(sorted(ties)))
