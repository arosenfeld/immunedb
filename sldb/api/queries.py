import datetime
import json
import math
import re
from sqlalchemy import desc, inspect, distinct
from sqlalchemy.sql import func
from sqlalchemy.ext.declarative import DeclarativeMeta

import sldb.util.lookups as lookups
from sldb.common.models import *
from sldb.common.mutations import MutationType, Mutations

def _fields_to_dict(fields, row):
    d = {}
    for f in fields:
        d[f] = getattr(row, f)
    return d


def _subject_to_dict(subject):
    return {
        'study': {
            'id': subject.study.id,
            'name': subject.study.name
        },
        'identifier': subject.identifier
    }

def _sample_to_dict(sample):
    d = _fields_to_dict(['id', 'name', 'info', 'subset', 'tissue',
                         'disease', 'lab', 'experimenter', 'valid_cnt',
                         'no_result_cnt', 'functional_cnt'], sample)
    d['date'] = sample.date.strftime('%Y-%m-%d')
    d['subject'] = _subject_to_dict(sample.subject)
    return d


def _clone_to_dict(clone):
    d = _fields_to_dict(['id', 'cdr3_nt'], clone)
    d['group'] = {
        'subject': _subject_to_dict(clone.group.subject),
        'v_gene': clone.group.v_gene,
        'j_gene': clone.group.j_gene,
        'cdr3_aa': clone.group.cdr3_aa,
        'cdr3_num_nts': clone.group.cdr3_num_nts,
        'germline': clone.group.germline,
    }
    return d


def get_all_studies(session):
    result = {}
    for sample in session.query(Sample).order_by(Sample.date):
        if sample.study.id not in result:
            result[sample.study.id] = {
                'id': sample.study.id,
                'name': sample.study.name,
                'info': sample.study.info,
                'samples': []
            }
        result[sample.study.id]['samples'].append(_sample_to_dict(sample))

    return result

def get_all_clones(session, filters, order_field, order_dir, paging=None):
    """Gets a list of all clones"""
    res = []
    clone_q = session.query(Clone)

    if filters is not None:
        for key, value in filters.iteritems():
            value = str(value).strip()
            if len(value) > 0 and value is not None:
                if key == 'min_cdr3_num_nts':
                    clone_q = clone_q.filter(Clone.cdr3_num_nts >= int(value))
                elif key == 'max_cdr3_num_nts':
                    clone_q = clone_q.filter(Clone.cdr3_num_nts <= int(value))
                else:
                    clone_q = clone_q.filter(
                        getattr(Clone, key).like(value.replace('*', '%')))

    if order_field is not None:
        order_field = getattr(Clone, order_field)
    else:
        order_field = Clone.id

    if order_dir is None or order_dir == 'desc':
        order_field = order_field.desc()
    else:
        order_field = order_field.asc()
        
    clone_q = clone_q.order_by(order_field)

    if paging is not None:
        page, per_page = paging
        clone_q = clone_q.offset((page - 1) * per_page).limit(per_page)

    for c in clone_q:
        stats_comb = []
        for stat in session.query(CloneFrequency)\
                           .filter(CloneFrequency.clone_id == c.id,
                                   CloneFrequency.filter_type == 'clones_all')\
                           .order_by(desc(CloneFrequency.total_sequences)):
                stats_comb.append({
                    'sample': {
                        'id': stat.sample.id,
                        'name': stat.sample.name
                    },
                    'unique_sequences': stat.unique_sequences,
                    'total_sequences': stat.total_sequences
                })

        clone_json = _clone_to_dict(c)
        clone_json['stats'] = stats_comb
        res.append(clone_json)

    return res


def compare_clones(session, uids):
    """Compares sequences within clones by determining their mutations"""
    clones = {}
    clone_muts = {}
    for clone_id, sample_id in uids:
        clone = session.query(Clone).filter(Clone.id == clone_id).first()
        germline = clone.germline[:309] + clone.cdr3_nt + \
            clone.germline[309 + clone.cdr3_num_nts:]
        if clone_id not in clones:
            clone_muts[clone_id] = Mutations(germline, clone.cdr3_num_nts)
            clone_json = _clone_to_dict(clone)
            clone_json['germline'] = germline
            clones[clone_id] = {
                'clone': clone_json,
                'mutation_stats': {},
                'seqs': []
            }
        mutations = clone_muts[clone_id]

        start_ptrn = re.compile('[N\-]*')
        for s in session.query(Sequence)\
                        .filter(Sequence.sample_id == sample_id,
                                Sequence.clone_id == clone_id)\
                        .group_by(Sequence.sequence_replaced):
            read_start = start_ptrn.match(s.sequence)
            if read_start is None:
                read_start = 0
            else:
                read_start = read_start.span()[1]

            clones[clone_id]['seqs'].append({
                'seq_id': s.seq_id,
                'sample': {
                    'id': s.sample.id,
                    'name': s.sample.name
                },
                'junction_nt': s.junction_nt,
                'sequence': s.sequence_replaced,
                'read_start': read_start,
                'mutations': mutations.add_sequence(s.sequence),
            })

        region_stats, pos_stats = mutations.get_aggregate()
        clones[clone_id]['mutation_stats']['regions'] = region_stats
        clones[clone_id]['mutation_stats']['positions'] = pos_stats

    return clones


def get_clone_overlap(session, filter_type, samples, subject=None,
                      paging=None):
    """Gets a list of clones and the samples in `samples` which they appear"""
    res = []
    q = session.query(
        CloneFrequency,
        func.sum(CloneFrequency.total_sequences).label('total_sequences'),
        func.group_concat(CloneFrequency.sample_id)
        .label('samples'))\
        .filter(CloneFrequency.filter_type == filter_type)

    if samples is not None:
        q = q.filter(CloneFrequency.sample_id.in_(samples))
    elif subject is not None:
        samples = map(
            lambda s: s.sample_id,
            session.query(distinct(Sequence.sample_id).label('sample_id'))
            .filter(Sequence.subject_id == subject).all())
        q = q.filter(CloneFrequency.sample_id.in_(samples))

    q = q.group_by(CloneFrequency.clone_id)\
        .order_by(desc(func.sum(CloneFrequency.total_sequences)))

    if paging is not None:
        page, per_page = paging
        # num_pages = math.ceil(len(q.all()) / per_page)
        q = q.offset((page - 1) * per_page).limit(per_page)

    for r in q:
        samples = ','.join(map(str, sorted(map(int, r.samples.split(',')))))
        freq = r.CloneFrequency
        unique_seqs = session.query(
            func.count(
                distinct(Sequence.sequence_replaced)).label('us')).filter(
            Sequence.clone == freq.clone).first().us
        res.append({
            'samples': samples,
            'unique_sequences': unique_seqs,
            'total_sequences': int(r.total_sequences),
            'clone': _clone_to_dict(freq.clone)
        })

    if paging:
        return res
    return res


def get_v_usage(session, filter_type, samples):
    """Gets the V-Gene usage percentages for samples"""
    data = {}
    headers = []
    for s in session.query(SampleStats)\
            .filter(SampleStats.filter_type == filter_type,
                    SampleStats.sample_id.in_(samples)):
        dist = json.loads(s.v_call_dist)
        data[s.sample.name] = {}
        total = 0
        for v in dist:
            total += v[1]

        for v in dist:
            name, occ = v
            if name in lookups.v_gene_ties:
                name = lookups.v_gene_ties[name]
            else:
                name = name.replace('/', '|').split('|')[0]
            if name not in headers:
                headers.append(name)

            data[s.sample.name][name] = round(100 * occ / float(total), 2)

    return data, headers


def get_all_subjects(session, paging):
    q = session.query(Subject)

    if paging is not None:
        page, per_page = paging
        q = q.offset((page - 1) * per_page).limit(per_page)

    subjects = []
    for subject in q:
        subjects.append({
            'id': subject.id,
            'identifier': subject.identifier,
            'study': {
                'id': subject.study.id,
                'name': subject.study.name
            },
            'total_samples': session.query(func.count(Sample.id)).filter(
                Sample.subject == subject).scalar(),
            'unique_seqs': session.query(func.count(Sequence.seq_id)).filter(
                Sequence.sample.has(subject=subject)).scalar(),
            'total_clones': session.query(func.count(Clone.id)).filter(
                Clone.subject_id == subject.id).scalar()
        })

    return subjects


def get_subject(session, sid):
    s = session.query(Subject).filter(Subject.id == sid).first()
    samples = []

    for seq in session.query(
            distinct(Sequence.sample_id).label('sample_id'))\
            .filter(Sequence.subject_id == s.id):
        # TODO: This should be a join but that is really slow for some reason
        sample = session.query(Sample).filter(
            Sample.id == seq.sample_id).first()

        samples.append({
            'id': sample.id,
            'name': sample.name,
            'date': sample.date.strftime('%Y-%m-%d'),
            'valid_cnt': sample.valid_cnt,
            'no_result_cnt': sample.no_result_cnt,
            'functional_cnt': sample.functional_cnt,
        })

    subject = {
        'id': s.id,
        'identifier': s.identifier,
        'study': {
            'id': s.study.id,
            'name': s.study.name,
        },
        'samples': samples,
        'unique_seqs': session.query(
            func.count(Sequence.seq_id).label('unique_seqs'))
        .filter(Sequence.subject_id == s.id).first().unique_seqs,
    }

    return subject


def get_sequence(session, sample_id, seq_id):
    seq = session.query(Sequence).filter(Sequence.sample_id == sample_id,
                                         Sequence.seq_id == seq_id).first()
    ret = _model_to_dict(seq)
    ret['subject'] = {
        'study': {
            'id': seq.subject.study.id,
            'name': seq.subject.study.name
        },
        'identifier': seq.subject.identifier,
    }
    ret['sample'] = {
        'id': seq.sample.id,
        'name': seq.sample.name,
    }

    if seq.clone is None:
        ret['clone'] = None
    else:
        ret['clone'] = _clone_to_dict(seq.clone)

    return ret
