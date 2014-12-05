import hashlib


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


def hash(seq_id, sample_id, sequence):
    return hashlib.sha1('{}{}{}'.format(
        seq_id, sample_id, sequence)).hexdigest()
