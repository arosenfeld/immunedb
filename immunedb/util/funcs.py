from collections import Counter


def consensus(strings):
    """Gets the unweighted consensus from a list of strings

    :param list strings: A set of equal-length strings.

    :returns: A consensus string
    :rtype: str

    """
    chrs = [Counter(chars).most_common(1)[0][0] for chars in zip(*strings)]
    return ''.join(chrs)


def get_regions(insertions):
    regions = [78, 36, 51, 30, 114]
    if insertions is not None and len(insertions) > 0:
        for pos, size in insertions:
            offset = 0
            for i, region_start in enumerate(regions):
                offset += region_start
                if pos < offset:
                    regions[i] += size
                    break

    return regions


def get_pos_region(regions, cdr3_len, pos):
    cdr3_start = sum(regions)
    j_start = cdr3_start + cdr3_len
    if pos >= j_start:
        return 'FR4'
    elif pos >= cdr3_start:
        return 'CDR3'

    total = 0
    for i, length in enumerate(regions):
        total += length
        if pos < total:
            rtype = 'FW' if i % 2 == 0 else 'CDR'
            rnum = (i // 2) + 1
            return '{}{}'.format(rtype, rnum)


def ord_to_quality(quality):
    if quality is None:
        return None
    return ''.join(map(lambda q: ' ' if q is None else chr(q + 33), quality))


def periodic_commit(session, query, interval=10):
    for i, r in enumerate(query):
        if i > 0 and i % interval == 0:
            session.commit()
        yield r
    session.commit()


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


def format_ties(ties):
    if ties is None:
        return None

    formatted = []
    for t in ties:
        prefix = t.prefix
        for e in t.name.split('|'):
            formatted.append(e.replace(prefix, '').split('*', 1)[0])
    return '{}{}'.format(prefix, '|'.join(sorted(set(formatted))))
