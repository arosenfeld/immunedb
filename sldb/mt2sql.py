import argparse
import sys
from os.path import basename

from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine
from sqlalchemy.engine.reflection import Inspector

from sldb.models import *
from sldb.stats import Stats


_mt_headers = ['order', 'seqID', 'functional', 'in-frame', 'stop',
               'mutation_invariate', 'v_match', 'v_length', 'j_match',
               'j_length', 'v_call', 'j_call', 'v_gap_length', 'j_gap_length',
               'juncton_gap_length', 'junction_nt', 'junction_aa',
               'gap_method', 'subject', 'subset', 'tissue', 'disease', 'date',
               'lab', 'experimenter', 'copy_number_close', 'collapse_to_close',
               'copy_number_iden', 'collapse_to_iden', 'sequence', 'germline',
               'cloneID', 'collapsedCloneID', 'withSinglecloneID']


_mt_field_map = {
    'seqID': 'seq_id',
    'in-frame': 'in_frame',

    'juncton_gap_length': None,
    'junction_aa': None,

    'germline': None,
    'cloneID': None,
    'collapsedCloneID': None,
    'withSinglecloneID': None
}


def _lookup_model_attrib(index):
    if _mt_headers[index] in _mt_field_map:
        return _mt_field_map[_mt_headers[index]]
    return _mt_headers[index]


def cast_to_type(table_type, field_name, value):
    ftype = getattr(Sequence.__table__.c, field_name).type
    if isinstance(ftype, Boolean):
        return value.upper() == 'T'
    elif isinstance(ftype, Integer):
        return int(value)
    elif isinstance(ftype, Date):
        month, day, year = map(int, value.split('/'))
        return datetime.datetime(year, month, day)
    return str(value)


def _get_or_create(session, model, **kwargs):
    ''' Gets or creates a record based on some kwargs search parameters '''
    instance = session.query(model).filter_by(**kwargs).first()
    if instance:
        return instance, False
    else:
        instance = model(**kwargs)
        session.add(instance)
        return instance, True


def _create_record(l):
    ''' Creates a new sequences record that has proper information '''
    record = Sequence()
    sp = map(lambda e: e.strip(), l.split('\t'))
    for i, val in enumerate(sp):
        if i < len(_mt_headers):
            model_attrib = _lookup_model_attrib(i)
            if model_attrib is not None and \
                    val not in ('unknown', 'NaN', ''):
                model_val = cast_to_type(Sequence, model_attrib, val)
                setattr(record, model_attrib, model_val)
    return record


def _create_noresult_record(l):
    ''' Creates a record for sequences which could not be assigned a V or J '''
    record = NoResult()
    seq_id = l.split('\t')[1].strip()  # TODO: This shouldn't be hardcoded
    record.seq_id = seq_id
    return record


def _get_clone(session, l):
    ''' Gets or creates a clone from a master table line '''
    sp = map(lambda e: e.strip(), l.split('\t'))
    clone = sp[_mt_headers.index('cloneID')]
    germline = sp[_mt_headers.index('germline')]
    if len(clone) == 0:
        return None, germline, None, None, None
    v, j, num_nts, cdr3, size, cn = clone.split('_')
    size = int(size)
    cn = int(cn)

    clone, new = _get_or_create(session, Clone, v_gene=v, j_gene=j, cdr3=cdr3,
                                cdr3_num_nts=num_nts)
    if new:
        clone.germline = germline
    else:
        if clone.germline != germline:
            print ('[ERROR] Germline mismatch with v={}, j={}, cdr3={},'
                   ' cdr3_len={}'
                   '\n\tIn DB    : {}'
                   '\n\tInserting: {}').format(
                       v, j, cdr3, num_nts, clone.germline, germline)
            sys.exit(1)
    return clone, germline, size, cn, new


def _n_to_germline(germ, seq, seq_id):
    if len(germ) != len(seq):
        print 'Germline and sequence different lengths'
        print 'SeqID', seq_id
        print germ
        print seq
        raise

    for i in range(0, len(seq)):
        if seq[i].upper() == 'N':
            seq = seq[:i] + germ[i] + seq[i+1:]
    return seq


def _add_mt(session, path, study_name, sample_name, sample_date, interval,
            force):
    ''' Processes the entire master table '''
    m, d, y = map(int, sample_date.split('-'))
    sample_date = '{}-{}-{}'.format(y, m, d)
    study, new = _get_or_create(session, Study, name=study_name)
    if new:
        print 'Created new study "{}"'.format(study.name)
    else:
        print 'Study "{}" already exists'.format(study.name)
    session.commit()

    print 'Adding {} sample for {} study'.format(sample_name, study_name)
    sample, new = _get_or_create(session, Sample,
                                 name=sample_name,
                                 study=study)

    if new:
        print '\tCreated new sample "{}"'.format(sample.name)
        sample.date = sample_date
        session.commit()
    else:
        if not force:
            print ('\tSample "{}" for study already exists.  '
                   'Skipping.').format(sample.name)
            return None
        else:
            print '\tForcing sample "{}"'.format(sample.name)
            sample = Sample(name=sample_name, study=study)
            sample.date = sample_date
            session.add(sample)
    session.commit()

    stats = Stats(session, sample.id)
    with open(path, 'rU') as fh:
        for i, l in enumerate(fh):
            pass
        total = i
        fh.seek(0)
        fh.readline()

        for i, l in enumerate(fh):
            if 'noresult' in l:
                record = _create_noresult_record(l)
            else:
                record = _create_record(l)

            record.sample = sample
            session.add(record)

            if 'noresult' not in l:
                clone, germline, size, cn, new = _get_clone(session, l)
                record.sequence_replaced = _n_to_germline(
                    germline,
                    record.sequence,
                    record.seq_id)
                if clone is not None:
                    record.clone = clone
                    record.clone_size = size
                    record.clone_copy_number = cn
                size = None
                cn = None
                stats.process_sequence(record, size, cn)
            else:
                stats.base_cnts['no_result_cnt'] += 1

            if i > 0 and i % interval == 0:
                session.commit()
                print '\t\tCommitted {} / {} ({}%)'.format(
                    i,
                    total,
                    round(100 * i / float(total), 2))
    session.commit()

    stats.add_and_commit()

    print 'Completed master table'
    return sample


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse master-table into \
    database.')
    parser.add_argument('host', help='mySQL host')
    parser.add_argument('db', help='mySQL database')
    parser.add_argument('user', help='mySQL user')
    parser.add_argument('pw', help='mySQL password')
    parser.add_argument('-c', type=int, default=1000, help='Number of'
                        ' sequences to parse between commits')
    parser.add_argument('-f', action='store_true', default=False, help='Forces'
                        ' insertion if sample name already exists')
    parser.add_argument('mt_dir', help='Master table directory')
    args = parser.parse_args()

    engine = create_engine(('mysql://{}:{}@{}/'
                            '{}?charset=utf8&use_unicode=0').format(
                                args.user, args.pw, args.host, args.db))

    Base.metadata.create_all(engine)
    Base.metadata.bind = engine
    session = sessionmaker(bind=engine)()

    for l in sys.stdin:
        study_name, sample_name, date, _ = map(lambda s: s.strip(),
                                               l.split('|'))
        path = '{}/{}/{}/{}'.format(
            args.mt_dir,  study_name, date, sample_name)
        print study_name, sample_name, date, path
        sample = _add_mt(
            session=session,
            path=path + '/master_table',
            study_name=study_name,
            sample_name=sample_name,
            sample_date=date,
            interval=args.c,
            force=args.f)
