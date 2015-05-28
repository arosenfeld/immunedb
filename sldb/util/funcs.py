import re

from sldb.common.models import Sequence


def trace_seq_collapses(session, seq):
    sample_col = session.query(
        Sequence.collapse_to_subject_seq_id,
        Sequence.collapse_to_subject_sample_id,
        Sequence.copy_number_in_sample
    ).filter(
        Sequence.seq_id == seq.collapse_to_sample_seq_id,
        Sequence.sample_id == seq.sample_id
    ).first()

    if sample_col is None:
        return None

    subject_col = session.query(
        Sequence.collapse_to_clone_seq_id,
        Sequence.collapse_to_clone_sample_id,
        Sequence.copy_number_in_subject
    ).filter(
        Sequence.seq_id == sample_col.collapse_to_subject_seq_id,
        Sequence.sample_id == sample_col.collapse_to_subject_sample_id,
    ).first()

    if subject_col is None:
        return None

    ret = {
        'sample': {
            'seq_id': seq.collapse_to_sample_seq_id,
            'sample_id': seq.sample_id,
            'copy_number': sample_col.copy_number_in_sample
        },
        'subject': {
            'seq_id': sample_col.collapse_to_subject_seq_id,
            'sample_id': sample_col.collapse_to_subject_sample_id,
            'copy_number': subject_col.copy_number_in_subject
        },
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


def page_query(q, per=10000, updating=False):
    """Pages a query in case it's too large for memory.  Extreme caution needs
    to be used if `updating` is True as rows may be returned multiple times
    unless the updates exclude the rows from further paging"""
    offset = 0
    while True:
        r = False
        q = q.limit(per)
        if not updating:
            q = q.offset(offset)
        for elem in q:
            r = True
            yield elem
        offset += per
        if not r:
            break


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
