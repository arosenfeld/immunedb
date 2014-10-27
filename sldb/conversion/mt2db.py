import argparse
import sys
import datetime
from os.path import basename
from collections import Counter

from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine, func
from sqlalchemy.engine.reflection import Inspector

import sldb.util.lookups as lookups
from sldb.common.models import *
from sldb.conversion.stats import Stats


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

    'germline': None,
    'cloneID': None,
    'collapsedCloneID': None,
    'withSinglecloneID': None
}

_cached_clones = {}

def _field_value(fields, key):
    return fields[_mt_headers.index(key)]


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


def _similarity(seq1, seq2):
    assert len(seq1) == len(seq2)
    return sum(ch1 == ch2 for ch1, ch2 in zip(seq1, seq2)) / float(len(seq1))


def _similar_to_all(seqs, check, min_similarity):
    for seq in seqs:
        if _similarity(check.junction_aa, seq.junction_aa) < min_similarity:
            return False
    return True


def _get_clone(session, seq, germline, min_similarity):
    ''' Gets or creates a clone from a master table line '''
    # Check if this exact CDR3 has already been assigned a clone
    key = (seq.v_call, seq.j_call, seq.junction_aa)
    if key in _cached_clones:
        c = _cached_clones[key]
    else:
        c = session.query(Clone).filter(
            Clone.v_gene == seq.v_call,
            Clone.j_gene == seq.j_call,
            Clone.cdr3_aa == seq.junction_aa).first()
        if c == None:
            # Otherwise, fuzzy match the CDR3
            for c in session.query(Clone).filter(
                Clone.v_gene==seq.v_call,
                Clone.j_gene==seq.j_call,
                Clone.cdr3_num_nts==len(seq.junction_nt)):

                if _similar_to_all(session.query(Sequence).filter(Sequence.clone == c),
                                   seq,
                                   min_similarity):
                    return c

            if seq.clone is None:
                c = Clone(
                    v_gene=seq.v_call,
                    j_gene=seq.j_call,
                    cdr3_nt=seq.junction_nt,
                    cdr3_num_nts=len(seq.junction_nt),
                    germline=germline)
    _cached_clones[key] = c
    return c


def _n_to_germline(seq, germline):
    if len(germline) != len(seq.sequence):
        print 'Germline and sequence different lengths'
        print 'SeqID', seq.seq_id
        print 'Germline: ', germline
        print 'Sequence: ', seq.sequence
        return seq

    filled_seq = ''
    for i, c in enumerate(seq.sequence):
        if c.upper() == 'N':
            filled_seq += germline[i].upper()
        else:
            filled_seq += c 
    return filled_seq


def _add_mt(session, path, study_name, sample_name, sample_date, interval,
            force, min_similarity):
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
    order_to_seq = {}
    with open(path, 'rU') as fh:
        fh.readline()
        for i, l in enumerate(fh):
            fields = l.split('\t')
            order = int(_field_value(fields, 'order'))
            seq = _field_value(fields, 'seqID')
            order_to_seq[order] = seq
        total = i
        fh.seek(0)
        fh.readline()

        for i, l in enumerate(fh):
            if 'noresult' in l:
                record = _create_noresult_record(l)
                stats.base_cnts['no_result_cnt'] += 1
                record.sample = sample
                session.add(record)
            else:
                record = _create_record(l)
                germline = _field_value(l.split('\t'), 'germline')
                record.sequence_replaced = _n_to_germline(record, germline)
                record.sample = sample
                if record.copy_number_iden > 0:
                    clone = _get_clone(session, record, germline,
                                       min_similarity)
                    '''
                    if len(clone.germline) != len(germline):
                        print ('Assigned sequence to clone with non-matching '
                               'germline length')
                        print 'Sample:', record.sample.name
                        print 'Seq ID:', record.seq_id
                        print 'Clone :', clone.germline
                        print 'Seq   :', germline
                    '''

                    record.clone = clone
                    stats.process_sequence(record)
                else:
                    record = DuplicateSequence(
                        sample=sample,
                        identity_seq_id=order_to_seq[record.collapse_to_iden],
                        seq_id=record.seq_id)

                session.add(record)

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


def _consensus(strings):
    cons = []
    for chars in zip(*strings):
        cons.append(Counter(chars).most_common(1)[0][0])

    return ''.join(cons)


def _process_all_clones(session):
    print 'Generating clone consensuses'
    for clone in session.query(Clone):
        clone.cdr3_nt = _consensus(
            map(lambda e: e.junction_nt,
            session.query(Sequence).filter(Sequence.clone_id==clone.id)))

        clone.cdr3_aa = lookups.aas_from_nts(
            clone.cdr3_nt[0:len(clone.cdr3_nt) - len(clone.cdr3_nt) % 3])
    session.commit()

def run_mt2db():
    parser = argparse.ArgumentParser(description='Parse master-table into \
    database.')
    parser.add_argument('host', help='mySQL host')
    parser.add_argument('db', help='mySQL database')
    parser.add_argument('user', help='mySQL user')
    parser.add_argument('pw', help='mySQL password')
    parser.add_argument('-s', type=int, default=85, help='Minimum similarity '
                        'between clone sequences.')
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
            force=args.f,
            min_similarity=args.s)

    _process_all_clones(session)
