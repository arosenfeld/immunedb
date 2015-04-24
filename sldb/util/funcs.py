import re

from sldb.common.models import Sequence


def collapse_seqs(session, seqs, copy_field, collapse_copy_field,
                  collapse_seq_id_field, collapse_sample_id_field=None):
    new_cns = {(s.sample_id, s.seq_id): getattr(s, copy_field) for s in seqs}

    for i, seq1 in enumerate(seqs):
        seq1_key = (seq1.sample_id, seq1.seq_id)
        if new_cns[seq1_key] == 0:
            continue
        pattern = seq_to_regex(seq1.sequence)
        for j, seq2 in enumerate(seqs[i+1:]):
            seq2_key = (seq2.sample_id, seq2.seq_id)
            if (new_cns[seq2_key] > 0 
                    and pattern.match(seq2.sequence) is not None):
                new_cns[seq1_key] += getattr(seq2, copy_field)
                new_cns[seq2_key] = 0 
                update_dict = {collapse_seq_id_field: seq1.seq_id}
                if collapse_sample_id_field is not None:
                    update_dict[collapse_sample_id_field] = seq1.sample_id
                session.query(Sequence).filter(
                    Sequence.sample_id == seq2.sample_id,
                    Sequence.seq_id == seq2.seq_id
                ).update(update_dict)

    for (sample_id, seq_id), cn in new_cns.iteritems():
        session.query(Sequence).filter(
            Sequence.sample_id == sample_id,
            Sequence.seq_id == seq_id
        ).update({
            collapse_copy_field: cn
        })


def seq_to_regex(seq):
    return re.compile(''.join(
        map(lambda c: '[{}N]'.format(c) if c != 'N' else '[ATCGN]', seq)
    ))


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
