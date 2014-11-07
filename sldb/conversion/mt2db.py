import argparse
import sys
import csv
import datetime
from os.path import basename
from collections import Counter

from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine, func
from sqlalchemy.engine.reflection import Inspector

import sldb.util.lookups as lookups
from sldb.common.models import *
from sldb.conversion.stats import Stats


_mt_remap_headers = {
    'seqID': 'seq_id',
    'in-frame': 'in_frame',
}


_mt_nomap = ['junction_gap_length', 'juncton_gap_length', 'subject',
             'germline', 'cloneID', 'collapsedCloneID', 'withSinglecloneID',
             'mutation_invariate']


_cached_clones = {}


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


def _create_sequence_record(row):
    ''' Creates a new sequences record that has proper information '''
    record = Sequence()
    for key, value in row.iteritems():
        if key not in _mt_nomap and value not in ('unknown', 'NaN', ''):
            setattr(record, key, cast_to_type(Sequence, key, value))

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
    key = (seq.v_call, seq.j_call, len(seq.junction_nt), seq.subject_id,
           seq.junction_aa)
    if key in _cached_clones:
        return _cached_clones[key]

    # Otherwise, fuzzy match the CDR3
    c = None
    for clone in session.query(Clone).filter(
            Clone.v_gene == seq.v_call,
            Clone.j_gene == seq.j_call,
            Clone.cdr3_num_nts == len(seq.junction_nt),
            Clone.subject == seq.subject):

        if _similar_to_all(session.query(Sequence).filter(
                Sequence.clone == clone), seq, min_similarity):
            c = clone
            break

    if c is None:
        c = Clone(
            v_gene=seq.v_call,
            j_gene=seq.j_call,
            cdr3_num_nts=len(seq.junction_nt),
            germline=germline,
            subject=seq.subject)
        session.add(c)

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


def _remap_headers(row):
    if None in row:
        del row[None]

    for remap_from, remap_to in _mt_remap_headers.iteritems():
        if remap_from in row:
            row[remap_to] = row.pop(remap_from)


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
    with open(path, 'rU') as fh:
        order_to_seq = {}
        for row in csv.DictReader(fh, delimiter='\t'):
            _remap_headers(row)
            order_to_seq[int(row['order'])] = row['seq_id']
        total = len(order_to_seq)
        fh.seek(0)

        for i, row in enumerate(csv.DictReader(fh, delimiter='\t')):
            _remap_headers(row)
            if 'noresult' in row.values():
                stats.base_cnts['no_result_cnt'] += 1
                session.add(NoResult(seq_id=row['seq_id'], sample=sample))
            else:
                record = _create_sequence_record(row)
                if record.copy_number_iden > 0:
                    record.subject, _ = _get_or_create(
                        session, Subject, study_id=study.id,
                        identifier=row['subject'])
                    record.sequence_replaced = _n_to_germline(record,
                                                              row['germline'])
                    record.sample = sample

                    if record.junction_nt is not None and \
                        record.junction_aa is not None and \
                        '*' not in record.junction_aa and \
                            record.copy_number_iden > 1:
                        clone = _get_clone(session, record, row['germline'],
                                           min_similarity)
                        if len(clone.germline) != len(row['germline']):
                            print '***Mismatched germline length***'
                            print 'Sample      :', sample.name
                            print 'Sequence    :', record.seq_id
                            print 'Length in DB:', len(clone.germline)
                            print 'Length in MT:', len(row['germline'])
                            clone = None
                        record.clone = clone

                    session.add(record)
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
            min_similarity=args.s / 100.0)
