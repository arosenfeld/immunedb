import datetime
import json
import math
import re
import numpy as np

from sqlalchemy import desc, inspect, distinct
from sqlalchemy.sql import func
from sqlalchemy.ext.declarative import DeclarativeMeta

import sldb.util.lookups as lookups
from sldb.common.models import *
from sldb.identification.vdj_sequence import VDJSequence
from sldb.common.mutations import MutationType, Mutations


_clone_filters = {
    'clones_all': lambda q: q,
    'clones_functional': lambda q: q.filter(
        SequenceMapping.functional == 1),
    'clones_nonfunctional': lambda q: q.filter(
        SequenceMapping.functional == 0),
}


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
        'identifier': subject.identifier,
        'id': subject.id
    }


def _sample_to_dict(sample):
    d = _fields_to_dict(['id', 'name', 'info', 'subset', 'tissue',
                         'disease', 'lab', 'experimenter'], sample)
    d['date'] = sample.date.strftime('%Y-%m-%d')
    d['subject'] = _subject_to_dict(sample.subject)
    return d


def _clone_to_dict(clone):
    d = _fields_to_dict(['id', 'cdr3_nt'], clone)
    d['group'] = {
        'id': clone.group.id,
        'v_gene': clone.group.v_gene,
        'j_gene': clone.group.j_gene,
        'cdr3_aa': clone.group.cdr3_aa,
        'cdr3_num_nts': clone.group.cdr3_num_nts,
        'subject': _subject_to_dict(clone.group.subject),
    }
    d['germline'] = clone.group.germline[0:VDJSequence.CDR3_OFFSET] + \
        clone.cdr3_nt + clone.group.germline[VDJSequence.CDR3_OFFSET +
                                             clone.group.cdr3_num_nts:]
    return d


def get_all_studies(session):
    result = {}
    for sample in session.query(Sample).order_by(Sample.date):
        if session.query(SequenceMapping).filter(
                SequenceMapping.sample == sample).first() is not None:
            status = 'reads'
        elif session.query(NoResult).filter(
                NoResult.sample == sample).first() is not None:
            status = 'noreads'
        else:
            status = 'unprocessed'

        if status in ('reads', 'noreads'):
            if sample.study.id not in result:
                result[sample.study.id] = {
                    'id': sample.study.id,
                    'name': sample.study.name,
                    'info': sample.study.info,
                    'samples': []
                }
            sample_dict = _sample_to_dict(sample)
            stats = session.query(SampleStats.sequence_cnt,
                                  SampleStats.in_frame_cnt,
                                  SampleStats.stop_cnt,
                                  SampleStats.functional_cnt,
                                  SampleStats.no_result_cnt).filter(
                SampleStats.sample_id == sample.id,
                SampleStats.outliers == 0,
                SampleStats.filter_type == 'all').first()
            if stats is not None:
                sample_dict['status'] = status
                sample_dict['sequence_cnt'] = stats.sequence_cnt
                sample_dict['in_frame_cnt'] = stats.in_frame_cnt
                sample_dict['stop_cnt'] = stats.stop_cnt
                sample_dict['functional_cnt'] = stats.functional_cnt
                sample_dict['no_result_cnt'] = stats.no_result_cnt
            result[sample.study.id]['samples'].append(sample_dict)

    return result


def get_all_clones(session, filters, order_field, order_dir, paging=None):
    """Gets a list of all clones"""
    res = []
    clone_q = session.query(Clone)

    if filters is not None:
        for key, value in filters.iteritems():
            if value is None:
                continue
            value = str(value).strip()
            if len(value) > 0 and value is not None:
                if key == 'min_cdr3_num_nts':
                    clone_q = clone_q.filter(Clone.cdr3_num_nts >= int(value))
                elif key == 'max_cdr3_num_nts':
                    clone_q = clone_q.filter(Clone.cdr3_num_nts <= int(value))
                elif key == 'id':
                    clone_q = clone_q.filter(Clone.id == int(value))
                elif key == 'group_id':
                    clone_q = clone_q.filter(Clone.group_id == int(value))
                else:
                    if hasattr(Clone, key):
                        c = Clone
                    else:
                        c = Clone.group
                    clone_q = clone_q.filter(
                        getattr(c, key).like(value.replace('*', '%')))

    if paging is not None:
        page, per_page = paging
        clone_q = clone_q.offset((page - 1) * per_page).limit(per_page)

    for c in clone_q:
        stats_comb = []
        for stat in session.query(
                SequenceMapping,
                func.count(SequenceMapping.identity_seq_id).label('unique'),
                func.sum(SequenceMapping.copy_number).label('total'))\
            .filter(SequenceMapping.clone_id == c.id)\
            .order_by(desc(func.count(
                SequenceMapping.identity_seq_id).label('unique')))\
                .group_by(SequenceMapping.sample_id):

                stats_comb.append({
                    'sample': {
                        'id': stat.SequenceMapping.sample.id,
                        'name': stat.SequenceMapping.sample.name
                    },
                    'unique_sequences': int(stat.unique),
                    'total_sequences': int(stat.total)
                })
        clone_dict = _clone_to_dict(c)
        clone_dict['stats'] = stats_comb
        res.append(clone_dict)

    return res


def compare_clones(session, uids):
    """Compares sequences within clones by determining their mutations"""
    clones = {}
    clone_muts = {}
    for clone_id, sample_ids in uids.iteritems():
        clone = session.query(Clone).filter(Clone.id == clone_id).first()
        germline = clone.group.germline
        if clone_id not in clones:
            clone_muts[clone_id] = Mutations(germline, clone.cdr3_nt)
            clones[clone_id] = {
                'clone': _clone_to_dict(clone),
                'mutation_stats': {},
                'seqs': []
            }
        mutations = clone_muts[clone_id]

        start_ptrn = re.compile('[N\-]*')

        q = session.query(SequenceMapping)\
            .join(Sequence)\
            .filter(SequenceMapping.clone_id == clone_id)
        if None not in sample_ids:
            q = q.filter(SequenceMapping.sample_id.in_(sample_ids))
        q = q.group_by(SequenceMapping.identity_seq_id)

        for mapping in q:
            read_start = start_ptrn.match(mapping.sequence)
            if read_start is None:
                read_start = 0
            else:
                read_start = read_start.span()[1]

            cn = session.query(
                func.sum(SequenceMapping.copy_number))\
                .filter(
                    SequenceMapping.identity_seq_id \
                        == mapping.identity_seq_id).scalar()
            muts = mutations.add_sequence(
                mapping.identity_seq.sequence_replaced)
            clones[clone_id]['seqs'].append({
                'seq_id': mapping.seq_id,
                'sample': {
                    'id': mapping.sample.id,
                    'name': mapping.sample.name,
                },
                'junction_nt': mapping.identity_seq.junction_nt,
                'sequence': mapping.identity_seq.sequence_replaced,
                'read_start': read_start,
                'copy_number': int(cn),
                'mutations': muts,
                'v_extent': mapping.v_length + \
                    mapping.num_gaps + mapping.pad_length,
                'j_length': mapping.j_length,
            })

        region_stats, pos_stats = mutations.get_aggregate()
        clones[clone_id]['mutation_stats']['regions'] = region_stats
        clones[clone_id]['mutation_stats']['positions'] = pos_stats

    return clones


def get_clone_tree(session, clone_id):
    return session.query(Clone.tree).filter(Clone.id == clone_id).first()


def get_clone_overlap(session, filter_type, ctype, limit,
                      paging=None):
    """Gets a list of clones and the samples in `samples` which they appear"""
    fltr = _clone_filters[filter_type]
    res = []
    q = fltr(session.query(
        SequenceMapping,
        func.group_concat(
            distinct(SequenceMapping.sample_id)).label('samples'),
        func.count(SequenceMapping.seq_id).label('unique'),
        func.sum(SequenceMapping.copy_number).label('total'))
        .join(Clone))

    if ctype == 'samples':
        q = q.filter(SequenceMapping.sample_id.in_(limit))
    elif ctype == 'subject':
        q = q.join(Sample).filter(Sample.subject_id == limit)

    q = q.order_by(desc('total')).group_by(SequenceMapping.clone_id)

    if paging is not None:
        page, per_page = paging
        q = q.offset((page - 1) * per_page).limit(per_page)

    for clone in q:
        res.append({
            'unique_sequences': int(clone.unique),
            'total_sequences': int(clone.total),
            'clone': _clone_to_dict(clone.SequenceMapping.clone),
            'samples': map(str, clone.samples.split(',')),
        })

    if paging:
        return res
    return res


def get_clones_in_samples(session, samples):
    return map(lambda e: e.id,
               session.query(
                   distinct(SequenceMapping.clone_id).label('id')).filter(
                   SequenceMapping.sample_id.in_(samples)))


def get_clones_in_subject(session, subject_id):
    return map(lambda e: e.id, session.query(Clone).filter(
        Clone.subject_id == subject_id))


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
            'unique_seqs': session.query(
                func.count(SequenceMapping.seq_id))
            .filter(SequenceMapping.sample.has(subject=subject)).scalar(),
            'total_clones': session.query(func.count(Clone.id)).filter(
                Clone.subject_id == subject.id).scalar()
        })

    return subjects


def get_subject(session, sid):
    s = session.query(Subject).filter(Subject.id == sid).first()
    samples = []

    for sample in session.query(Sample).filter(Sample.subject_id == sid):
        stats = session.query(SampleStats).filter(
            SampleStats.filter_type == 'all',
            SampleStats.sample_id == sample.id,
            SampleStats.outliers == True).first()
        sample_dict = {
            'id': sample.id,
            'name': sample.name,
            'date': sample.date.strftime('%Y-%m-%d'),
        }
        if stats is not None:
            sample_dict['valid_cnt'] = stats.sequence_cnt
            sample_dict['no_result_cnt'] = stats.no_result_cnt
            sample_dict['functional_cnt'] = stats.functional_cnt
        samples.append(sample_dict)

    subject = {
        'id': s.id,
        'identifier': s.identifier,
        'study': {
            'id': s.study.id,
            'name': s.study.name,
        },
        'samples': samples,
    }

    return subject


def get_stats(session, samples, include_outliers):
    counts = {}
    stats = {}
    dist_fields = [
        'v_match_dist', 'v_length_dist', 'j_match_dist',
        'j_length_dist', 'v_call_dist', 'j_call_dist',
        'cdr3_length_dist', 'copy_number_dist']
    cnt_fields = ['sequence_cnt', 'in_frame_cnt', 'stop_cnt', 'functional_cnt',
                 'no_result_cnt']
    for stat in session.query(SampleStats).filter(
            SampleStats.sample_id.in_(samples),
            SampleStats.outliers == include_outliers):
        if stat.sample_id not in stats:
            stats[stat.sample_id] = {
                'sample': _sample_to_dict(stat.sample),
                'filters': {},
            }
        if stat.filter_type not in counts:
            counts[stat.filter_type] = 0

        flds = _fields_to_dict(dist_fields + cnt_fields, stat)
        stats[stat.sample_id]['filters'][stat.filter_type] = flds
        counts[stat.filter_type] += stat.sequence_cnt

    return {'counts': counts, 'stats': stats }


def get_sequence(session, sample_id, seq_id):
    seq = session.query(SequenceMapping)\
        .filter(SequenceMapping.sample_id == sample_id,
                SequenceMapping.seq_id == seq_id).first()
    if seq is None:
        seq = session.query(DuplicateSequence)\
            .filter(DuplicateSequence.seq_id == seq_id).first()

        seq = session.query(SequenceMapping)\
            .filter(SequenceMapping.sample_id == sample_id,
                    SequenceMapping.seq_id == seq.seq_id).first()

    ret = _fields_to_dict(['seq_id', 'identity_seq_id', 'alignment', 'v_match',
                           'j_match', 'v_length', 'j_length', 'in_frame',
                           'functional', 'stop', 'copy_number', 'sequence',
                           'pre_cdr3_length', 'pre_cdr3_match',
                           'post_cdr3_length', 'post_cdr3_match', 'pad_length',
                           'num_gaps', 'levenshtein_dist'],
                          seq)
    ret['sample'] = _sample_to_dict(seq.sample)

    ret['v_call'] = seq.identity_seq.v_call
    ret['j_call'] = seq.identity_seq.j_call
    ret['junction_nt'] = seq.identity_seq.junction_nt
    ret['junction_aa'] = seq.identity_seq.junction_aa
    ret['germline'] = seq.identity_seq.germline
    ret['v_extent'] = ret['v_length'] + ret['num_gaps'] + ret['pad_length']

    if seq.clone is None:
        ret['clone'] = None
    else:
        ret['clone'] = _clone_to_dict(seq.clone)

    muts = Mutations(seq.identity_seq.germline,
                     seq.identity_seq.junction_nt)
    ret['mutations'] = muts.add_sequence(seq.sequence)

    ret['duplicates'] = []
    ret['total_copy_number'] = ret['copy_number']
    for dup in session.query(SequenceMapping).filter(
            SequenceMapping.identity_seq_id == seq.identity_seq_id,
            SequenceMapping.sample.has(subject_id=seq.sample.subject_id),
            SequenceMapping.seq_id != seq_id)\
            .order_by(SequenceMapping.sample_id):
        ret['duplicates'].append({
            'seq_id': dup.seq_id,
            'sample': _sample_to_dict(dup.sample),
            'alignment': dup.alignment,
            'copy_number': dup.copy_number,
        })
        ret['total_copy_number'] += dup.copy_number

    return ret


def get_all_sequences(session, filters, order_field, order_dir, paging=None):
    """Gets a list of all clones"""
    def get_field(key):
        tbls = [SequenceMapping, Sequence, Subject, Clone]
        for t in tbls:
            if hasattr(t, key):
                return getattr(t, key)

    res = []
    query = session.query(SequenceMapping).join(Sequence).join(Sample)\
        .outerjoin(Clone)

    if filters is not None:
        for key, value in filters.iteritems():
            if value in [None, True, False]:
                continue
            value = str(value).strip()
            if len(value) > 0 and value is not None:
                if key == 'sample_id':
                    query = query.filter(SequenceMapping.sample_id ==
                                         int(value))
                elif key == 'in_frame':
                    query = query.filter(
                        SequenceMapping.in_frame == int(value))
                elif key == 'min_copy_number':
                    query = query.filter(
                        SequenceMapping.copy_number >= int(value))
                elif key == 'max_copy_number':
                    query = query.filter(
                        SequenceMapping.copy_number <= int(value))
                else:
                    query = query.filter(get_field(key).like(
                        value.replace('*', '%')))

    if filters is None or 'show_r1' not in filters or not filters['show_r1']:
        query = query.filter(SequenceMapping.alignment == 'pRESTO')
    if filters is None or 'show_indel' not in filters or not filters['show_indel']:
        query = query.filter(SequenceMapping.levenshtein_dist.is_(None))

    if paging is not None:
        page, per_page = paging
        query = query.offset((page - 1) * per_page).limit(per_page)

    for row in query:
        fields = _fields_to_dict(
            ['seq_id', 'alignment', 'v_match', 'j_match', 'v_length',
             'j_length', 'in_frame', 'functional', 'stop',
             'levenshtein_dist'], row)

        fields['copy_number'] = int(session.query(
            func.sum(SequenceMapping.copy_number)).filter(
                SequenceMapping.unique_id == row.unique_id).scalar())
        fields = dict(fields.items() + _fields_to_dict(
            ['v_call', 'j_call', 'junction_aa'], row.identity_seq).items())

        fields['sample'] = _sample_to_dict(row.sample)
        fields['cdr3_length'] = len(row.identity_seq.junction_nt)
        if row.clone is None:
            fields['clone'] = None
        else:
            fields['clone'] = _clone_to_dict(row.clone)
        res.append(fields)

    return res

_export_fields = [
    'seq_id',
    'identity_seq_id',
    ('subject_id', lambda seq: seq.sample.subject.id),
    ('subject_identifier', lambda seq: seq.sample.subject.identifier),
    ('subset', lambda seq: seq.sample.subset),
    ('tissue', lambda seq: seq.sample.tissue),
    ('disease', lambda seq: seq.sample.disease),
    ('lab', lambda seq: seq.sample.lab),
    ('experimenter', lambda seq: seq.sample.experimenter),
    ('date', lambda seq: seq.sample.date),
    
    'sample_id',
    ('sample_name', lambda seq: seq.sample.name),
    ('study_id', lambda seq: seq.sample.study.id),
    ('study_name', lambda seq: seq.sample.study.name),

    'alignment',
    'levenshtein_dist',

    'num_gaps',
    'pad_length',
    
    'v_match',
    'v_length',
    'j_match',
    'j_length',

    'pre_cdr3_match',
    'pre_cdr3_length',
    'post_cdr3_match',
    'post_cdr3_length',

    'in_frame',
    'functional',
    'stop',
    'copy_number',
    
    'sequence',
    ('sequence_filled', lambda seq: seq.identity_seq.sequence_replaced),
    ('germline', lambda seq: seq.identity_seq.germline),

    ('v_call', lambda seq: seq.identity_seq.v_call),
    ('j_call', lambda seq: seq.identity_seq.j_call),
    ('cdr3_nt', lambda seq: seq.identity_seq.junction_nt),
    ('cdr3_aa', lambda seq: seq.identity_seq.junction_aa),
    ('gap_method', lambda seq: seq.identity_seq.gap_method),

    'clone_id',
    ('clone_group_id', lambda seq: seq.clone.group.id),
    ('clone_cdr3_nt', lambda seq: seq.clone.cdr3_nt),
    ('clone_cdr3_aa', lambda seq: seq.clone.group.cdr3_aa),
    ('clone_json_tree', lambda seq: seq.clone.tree),
]


def _name_and_field(e):
    if type(e) == str:
        return e, lambda seq: getattr(seq, e)
    return e


def _fasta_header(no_key, keys):
    return '>{}{}{}\n'.format(
        no_key,
        '|' if len(no_key) > 0 and len(keys) > 0 else '',
        '|'.join(map(lambda (k, v): '{}={}'.format(k, v), keys)))


def export_seqs(session, eformat, rtype, rids, selected_fields,
                duplicates, noresults):
    seqs = session.query(SequenceMapping).filter(
        getattr(SequenceMapping, '{}_id'.format(rtype)).in_(rids))
    if eformat == 'clip':
        seqs = seqs.order_by(SequenceMapping.identity_seq_id)

    if eformat == 'tab':
        headers = []
        for field in _export_fields:
            n, f = _name_and_field(field)
            if n in selected_fields:
                headers.append(n)

        yield '{}\n'.format('\t'.join(headers))

    last_iden = None
    for seq in seqs:
        data = []
        for field in _export_fields:
            n, f = _name_and_field(field)
            if n in selected_fields:
                if 'clone' not in n or (seq.clone is not None):
                    data.append((n, f(seq)))
                else:
                    data.append((n, 'NA'))

        if eformat == 'tab':
            yield '{}\n'.format('\t'.join(map(lambda s: str(s[1]), 
                                              data)))
        else:
            if eformat == 'clip' and last_iden != seq.identity_seq_id:
                last_iden = seq.identity_seq_id
                yield _fasta_header(
                    '>Germline',
                    (('v_gene', seq.identity_seq.v_call),
                     ('j_gene', seq.identity_seq.j_call),
                     ('cdr3_aa', seq.identity_seq.junction_aa),
                     ('cdr3_len', seq.identity_seq.junction_num_nts)))

            yield _fasta_header(seq.seq_id, data)
            yield '{}\n'.format(seq.identity_seq.sequence_replaced \
                if eformat == 'fill' else seq.sequence)
