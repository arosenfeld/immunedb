from collections import Counter
import itertools

import dnautils


def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]


def yield_limit(qry, pk_attr, maxrq=5000):
    firstid = None
    while True:
        q = qry
        if firstid is not None:
            q = qry.filter(pk_attr > firstid)
        rec = None
        for rec in q.order_by(pk_attr).limit(maxrq):
            yield rec
        if rec is None:
            break
        firstid = pk_attr.__get__(rec, pk_attr) if rec else None


def bulk_add(session, objs, chunk_size=100, flush=True):
    for i in range(0, len(objs), chunk_size):
        session.bulk_insert_mappings(
            type(objs[0]),
            [o.__dict__ for o in objs[i:i + chunk_size]]
        )
        if flush:
            session.flush()


def flatten(iterable):
    return list(itertools.chain.from_iterable(iterable))


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
    return ''.join([' ' if q is None else chr(q + 33) for q in quality])


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


def format_ties(ties, strip_alleles=True):
    if ties is None:
        return None

    formatted = []
    for t in ties:
        prefix = t.prefix
        for e in t.name.split('|'):
            name = e.replace(prefix, '')
            if strip_alleles:
                name = name.split('*', 1)[0]
            formatted.append(name)
    return '{}{}'.format(prefix, '|'.join(sorted(set(formatted))))


def get_cigar(ref, qry):
    assert len(ref) == len(qry)
    cnt = 0
    current_op = None
    cigar = []
    for r, q in zip(ref, qry):
        if r == q == '-':
            continue
        elif dnautils.equal(r, q):
            op = '='
        elif r == '-':
            op = 'I'
        elif q == '-':
            op = 'D'
        else:
            op = 'X'
        if op != current_op:
            if current_op:
                cigar.append('{}{}'.format(cnt, current_op))
            current_op = op
            cnt = 1
        else:
            cnt += 1
    cigar.append('{}{}'.format(cnt, current_op))
    return ''.join(cigar)
