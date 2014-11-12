def get_or_create(session, model, **kwargs):
    ''' Gets or creates a record based on some kwargs search parameters '''
    instance = session.query(model).filter_by(**kwargs).first()
    if instance:
        return instance, False
    else:
        instance = model(**kwargs)
        session.add(instance)
        return instance, True
