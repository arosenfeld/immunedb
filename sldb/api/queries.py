import datetime
import json
import math
import numpy as np
import re

from sqlalchemy import desc, distinct, inspect, or_
from sqlalchemy.sql import func
from sqlalchemy.sql.expression import false, true
from sqlalchemy.ext.declarative import DeclarativeMeta

import sldb.util.lookups as lookups
import sldb.util.funcs as funcs
from sldb.common.models import *
from sldb.common.mutations import CloneMutations, threshold_mutations
from sldb.identification.v_genes import VGene


_clone_filters = {
    'clones_all': lambda q: q,
    'clones_functional': lambda q: q.filter(
        Clone.cdr3_num_nts % 3 == 0),
    'clones_nonfunctional': lambda q: q.filter(
        Clone.cdr3_num_nts % 3 != 0),
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
    d = _fields_to_dict(['id', 'cdr3_nt', 'v_gene', 'j_gene', 'cdr3_aa',
                         'cdr3_num_nts'], clone)
    d['subject'] = _subject_to_dict(clone.subject)
    d['germline'] = clone.consensus_germline

    return d


def _page_query(q, paging):
    if paging is None:
        return q
    page, per_page = paging
    return q.offset((page - 1) * per_page).limit(per_page)


def get_all_studies(session):
    result = {}
    for sample in session.query(Sample).order_by(Sample.date):
        if session.query(Sequence.seq_id).filter(
                Sequence.sample == sample).first() > 0:
            status = 'reads'
        elif session.query(NoResult.seq_id).filter(
                NoResult.sample == sample).first() > 0:
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
                SampleStats.outliers == true(),
                SampleStats.full_reads == false(),
                SampleStats.filter_type == 'all').first()
            if stats is not None:
                sample_dict['status'] = status
                sample_dict['sequence_cnt'] = stats.sequence_cnt
                sample_dict['in_frame_cnt'] = stats.in_frame_cnt
                sample_dict['stop_cnt'] = stats.stop_cnt
                sample_dict['functional_cnt'] = stats.functional_cnt
                sample_dict['no_result_cnt'] = stats.no_result_cnt
            else:
                sample_dict['status'] = 'processing'
            result[sample.study.id]['samples'].append(sample_dict)

    return result


def get_all_clones(session, filters, order_field, order_dir, paging=None):
    """Gets a list of all clones"""
    def get_field(key):
        tbls = [Clone, CloneStats]
        for t in tbls:
            if hasattr(t, key):
                return getattr(t, key)

    res = []
    clone_q = session.query(
        Clone, CloneStats.unique_cnt, CloneStats.total_cnt
    ).join(CloneStats).filter(
        CloneStats.sample_id.is_(None)
    )
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
                elif key == 'min_unique':
                    clone_q = clone_q.filter(
                        CloneStats.unique_cnt >= int(value))
                elif key == 'max_unique':
                    clone_q = clone_q.filter(
                        CloneStats.unique_cnt <= int(value))
                elif key == 'id':
                    clone_q = clone_q.filter(Clone.id == int(value))
                else:
                    clone_q = clone_q.filter(getattr(Clone, key).like(
                        value.replace('*', '%')))

    order_field = get_field(order_field)

    if order_dir == 'asc':
        clone_q = clone_q.order_by(order_field)
    else:
        clone_q = clone_q.order_by(desc(order_field))

    for c, unique_cnt, total_cnt in _page_query(clone_q, paging):
        stats_comb = []

        query = session.query(
            CloneStats.unique_cnt, CloneStats.total_cnt, Sample
        ).join(Sample).filter(
            CloneStats.clone_id == c.id,
            ~CloneStats.sample_id.is_(None)
        ).order_by(desc(CloneStats.unique_cnt))
        for stat in query:
            stats_comb.append({
                'sample': {
                    'id': stat.Sample.id,
                    'name': stat.Sample.name
                },
                'unique_sequences': int(stat.unique_cnt),
                'total_sequences': int(stat.total_cnt)
            })
        clone_dict = _clone_to_dict(c)
        clone_dict['unique_sequences'] = unique_cnt
        clone_dict['total_sequences'] = total_cnt
        clone_dict['stats'] = stats_comb
        res.append(clone_dict)

    return res


def get_clone(session, clone_id, sample_ids, thresholds=None):
    """Compares sequences within clones by determining their mutations"""

    if thresholds is None:
        thresholds = [
            ('percent', 100),
            ('percent', 80),
            ('percent', 50),
            ('percent', 20),
            ('percent', 0),
            ('seqs', 2),
            ('seqs', 5),
            ('seqs', 10),
            ('seqs', 25)
        ]

    result = {}
    clone = session.query(Clone).filter(Clone.id == clone_id).first()

    q = session.query(
        Sequence,
        Sequence.copy_number_in_subject
    ).filter(
        Sequence.clone_id == clone_id,
        Sequence.copy_number_in_subject > 0
    )
    if sample_ids is not None:
        q = q.filter(Sequence.sample_id.in_(sample_ids))

    result = {
        'clone': _clone_to_dict(clone),
        'quality': [],
        'seqs': []
    }

    start_ptrn = re.compile('[N\-]*')
    for seqr in q:
        seq = seqr.Sequence
        read_start = start_ptrn.match(seq.sequence).span()[1] or 0
        if seq.quality is not None:
            diff = len(seq.quality) - len(result['quality'])
            if diff > 0:
                result['quality'].extend([[] for _ in range(0, diff)])
            for i, b in enumerate(seq.quality):
                if b is not ' ':
                    result['quality'][i].append(ord(b) - 33)

        result['seqs'].append({
            'seq_id': seq.seq_id,
            'sample': {
                'id': seq.sample.id,
                'name': seq.sample.name,
            },
            'cdr3_nt': seq.cdr3_nt,
            'sequence': seq.sequence,
            'read_start': read_start,
            'copy_number_in_subject': int(seqr.copy_number_in_subject),
            'mutations': json.loads(seq.mutations_from_clone),
            'v_extent': seq.v_length + seq.num_gaps + seq.pad_length,
            'j_length': seq.j_length,
        })

    res_qual = []
    for i, quals in enumerate(result['quality']):
        if len(quals) > 0:
            res_qual.append((i, round(np.mean(quals), 2)))
    result['quality'] = res_qual

    all_mutations, total_seqs = get_clone_mutations(session, clone_id,
                                                    sample_ids)
    mut_dict = {
        'positions': all_mutations['positions'],
        'regions': {}
    }
    for threshold in thresholds:
        if threshold[0] == 'seqs':
            seq_min = threshold[1]
        else:
            seq_min = int(math.ceil(threshold[1] / 100.0 * total_seqs))
        tname = '_'.join(map(str, threshold))
        mut_dict['regions'][tname] = threshold_mutations(all_mutations,
                                                         seq_min)
    result['mutation_stats'] = mut_dict

    stats = session.query(
        CloneStats
    ).filter(
        CloneStats.clone_id == clone.id,
        ~CloneStats.sample_id.is_(None)
    ).order_by(desc(CloneStats.unique_cnt))
    result['samples'] = []
    for stat in stats:
        sample = _sample_to_dict(stat.sample)
        sample['unique'] = stat.unique_cnt
        sample['total'] = stat.total_cnt
        result['samples'].append(sample)

    return result


def get_selection_pressure(session, clone_id, sample_ids):
    query = session.query(
        CloneStats.sample_id, CloneStats.selection_pressure
    ).filter(CloneStats.clone_id == clone_id)
    if sample_ids is not None:
        query = query.filter(
            or_(
                CloneStats.sample_id.in_(sample_ids),
                CloneStats.sample_id.is_(None)
            )
        )

    pressure = []
    for row in query:
        if row.sample_id is None:
            name = 'All'
        else:
            name = session.query(Sample.name).filter(
                Sample.id == row.sample_id).first().name
        pressure.append({
            'sample': {
                'id': row.sample_id,
                'name': name
            },
            'pressure': json.loads(row.selection_pressure)
        })

    return pressure


def get_clone_mutations(session, clone_id, sample_ids=None):
    if sample_ids is None or type(sample_ids) is int or len(sample_ids) == 1:
        if sample_ids is None:
            sample_id = None
        elif type(sample_ids) is int:
            sample_id = sample_ids
        elif len(sample_ids) == 1:
            sample_id = sample_ids[0]

        clone_stats = session.query(
            CloneStats.mutations,
            CloneStats.unique_cnt
        ).filter(
            CloneStats.clone_id == clone_id,
            CloneStats.sample_id == sample_id
        ).first()
        all_mutations = json.loads(clone_stats.mutations)
        total_seqs = clone_stats.unique_cnt
    else:
        clone = session.query(Clone).filter(Clone.id == clone_id).first()
        all_mutations = CloneMutations(session, clone).calculate(
            limit_samples=sample_ids
        ).get_all()
        total_seqs = session.query(func.count(Sequence.seq_id)).filter(
            Sequence.clone_id == clone_id,
            Sequence.sample_id.in_(sample_ids)).scalar()

    return all_mutations, total_seqs


def get_clone_tree(session, clone_id):
    return session.query(Clone.tree).filter(Clone.id == clone_id).first()


def get_clone_overlap(session, filter_type, ctype, limit,
                      paging=None):
    """Gets a list of clones and the samples in `samples` which they appear"""
    fltr = _clone_filters[filter_type]
    res = []

    if ctype == 'samples':
        clones = session.query(
            Clone,
            func.sum(CloneStats.unique_cnt).label('unique_cnt'),
            func.sum(CloneStats.total_cnt).label('total_cnt')
        ).join(CloneStats).filter(
            CloneStats.sample_id.in_(limit),
        ).group_by(CloneStats.clone_id)
    elif ctype == 'subject':
        clones = session.query(
            Clone, CloneStats.unique_cnt.label('unique_cnt'),
            CloneStats.total_cnt.label('total_cnt')
        ).join(CloneStats).filter(
            CloneStats.sample_id.is_(None),
            Clone.subject_id == limit
        )

    clones = fltr(clones.order_by(desc('unique_cnt')))

    for clone in _page_query(clones, paging):
        selected_samples = []
        other_samples = []
        query = session.query(CloneStats).filter(
            CloneStats.clone_id == clone.Clone.id,
            ~CloneStats.sample_id.is_(None)
        ).order_by(
            desc(CloneStats.total_cnt)
        )
        for stat in query:
            data = {
                'id': stat.sample_id,
                'name': stat.sample.name,
                'unique_sequences': stat.unique_cnt,
                'total_sequences': stat.total_cnt
            }
            if ctype == 'subject' or stat.sample_id in limit:
                selected_samples.append(data)
            else:
                other_samples.append(data)

        res.append({
            'unique_sequences': int(clone.unique_cnt),
            'total_sequences': int(clone.total_cnt),
            'clone': _clone_to_dict(clone.Clone),
            'selected_samples': selected_samples,
            'other_samples': other_samples,
        })

    if paging:
        return res
    return res


def get_clones_in_samples(session, samples):
    return map(lambda e: e.id,
               session.query(
                   distinct(Sequence.clone_id).label('id')).filter(
                   Sequence.sample_id.in_(samples)))


def get_clones_in_subject(session, subject_id):
    return map(lambda e: e.id, session.query(Clone).filter(
        Clone.subject_id == subject_id))


def get_v_usage(session, samples, filter_type, include_outliers,
                include_partials, grouping, by_family):
    """Gets the V-Gene usage percentages for samples"""
    if by_family:
        def name_func(s):
            return s.split('*')[0].split('-', 1)[0].replace('IGHV', '')
    else:
        def name_func(s):
            return s.split('*')[0].replace('IGHV', '')
    data = {}
    totals = {}
    for s in session.query(SampleStats)\
            .filter(SampleStats.filter_type == filter_type,
                    SampleStats.outliers == include_outliers,
                    SampleStats.full_reads != include_partials,
                    SampleStats.sample_id.in_(samples)):
        dist = json.loads(s.v_gene_dist)
        if grouping == 'subject':
            group_key = s.sample.subject.identifier
        else:
            group_key = getattr(s.sample, grouping)

        if group_key is None:
            group_key = 'None'

        if group_key not in data:
            data[group_key] = {}

        for v in dist:
            name, occ = v
            name = '|'.join(sorted(set(map(name_func, name.split('|')))))

            if name not in data[group_key]:
                data[group_key][name] = 0
            data[group_key][name] += occ

    headers = []
    for group_key, names in data.iteritems():
        totals[group_key] = float(sum(names.values()))
        for name, value in names.iteritems():
            percent = round(100 * value / totals[group_key], 2)
            names[name] = percent
            if name not in headers and percent >= 1.0:
                headers.append(name)

    return data, sorted(headers), totals


def get_all_subjects(session, paging=None):
    subjects = []
    for subject in _page_query(session.query(Subject), paging):
        seqs = session.query(func.sum(SampleStats.sequence_cnt))\
            .filter(
                SampleStats.sample.has(subject=subject),
                SampleStats.filter_type == 'all',
                SampleStats.outliers == true(),
                SampleStats.full_reads == false()).scalar()
        info = {
            'id': subject.id,
            'identifier': subject.identifier,
            'study': {
                'id': subject.study.id,
                'name': subject.study.name
            },
        }
        if seqs is not None:
            info['total_samples'] = session.query(
                func.count(Sample.id)
            ).filter(
                Sample.subject == subject
            ).scalar()
            info['unique_seqs'] = int(seqs)
            info['total_clones'] = session.query(
                func.count(Clone.id)
            ).filter(
                Clone.subject_id == subject.id
            ).scalar()
        subjects.append(info)

    return subjects


def get_subject(session, sid):
    s = session.query(Subject).filter(Subject.id == sid).first()
    samples = []

    for sample in session.query(Sample).filter(Sample.subject_id == sid):
        stats = session.query(SampleStats).filter(
            SampleStats.filter_type == 'all',
            SampleStats.sample_id == sample.id,
            SampleStats.outliers == true(),
            SampleStats.full_reads == false()).first()
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


def get_stats(session, samples, include_outliers, include_partials, grouping):
    counts = {}
    stats = {}
    sample_info = {}
    dist_fields = [
        'v_match_dist', 'v_length_dist', 'v_identity_dist', 'j_match_dist',
        'j_length_dist', 'v_gene_dist', 'j_gene_dist', 'cdr3_length_dist',
        'copy_number_dist', 'quality_dist'
    ]
    cnt_fields = [
        'sequence_cnt', 'in_frame_cnt', 'stop_cnt', 'functional_cnt',
        'no_result_cnt'
    ]

    for stat in session.query(SampleStats).filter(
            SampleStats.sample_id.in_(samples),
            SampleStats.outliers == include_outliers,
            SampleStats.full_reads != include_partials):
        if grouping == 'subject':
            group_key = stat.sample.subject.identifier
        else:
            group_key = getattr(stat.sample, grouping)

        if group_key not in stats:
            stats[group_key] = {}

        if stat.sample.id not in sample_info:
            sample_info[stat.sample.id] = _sample_to_dict(stat.sample)

        if stat.filter_type not in counts:
            counts[stat.filter_type] = {'total': 0}
        if stat.sample.id not in counts[stat.filter_type]:
            counts[stat.filter_type][stat.sample.id] = 0
        counts[stat.filter_type][stat.sample.id] = stat.sequence_cnt
        counts[stat.filter_type]['total'] += stat.sequence_cnt

        fields = _fields_to_dict(dist_fields + cnt_fields, stat)
        if stat.filter_type not in stats[group_key]:
            stats[group_key][stat.filter_type] = {}

        for field, values in fields.iteritems():
            if field.endswith('cnt'):
                if stat.filter_type == 'all':
                    sample_info[stat.sample.id][field] = values
                continue
            if field not in stats[group_key][stat.filter_type]:
                stats[group_key][stat.filter_type][field] = {}

            for (x, freq) in json.loads(values):
                if x not in stats[group_key][stat.filter_type][field]:
                    stats[group_key][stat.filter_type][field][x] = 0
                stats[group_key][stat.filter_type][field][x] += freq

    for group, filter_dict in stats.iteritems():
        for filter_name, key_dict in filter_dict.iteritems():
            for key, vals in key_dict.iteritems():
                reduced = []
                for x in sorted(vals.keys()):
                    reduced.append((x, vals[x]))
                key_dict[key] = reduced

    return {'samples': sample_info, 'counts': counts, 'stats': stats}


def get_sequence(session, sample_id, seq_id):
    seq = session.query(Sequence)\
        .filter(Sequence.sample_id == sample_id,
                Sequence.seq_id == seq_id).first()
    if seq is None:
        dup_seq = session.query(DuplicateSequence.duplicate_seq_id)\
            .filter(DuplicateSequence.seq_id == seq_id).first()
        seq = session.query(Sequence)\
            .filter(Sequence.sample_id == sample_id,
                    Sequence.seq_id == dup_seq.duplicate_seq_id).first()

    ret = _fields_to_dict([
        'seq_id', 'partial', 'paired', 'v_gene', 'j_gene', 'cdr3_nt',
        'cdr3_aa', 'germline', 'v_match', 'j_match', 'v_length',
        'j_length', 'in_frame', 'functional', 'stop', 'copy_number',
        'sequence', 'pre_cdr3_length', 'pre_cdr3_match', 'post_cdr3_length',
        'post_cdr3_match', 'pad_length', 'num_gaps',
        'probable_indel_or_misalign', 'quality',
    ], seq)
    ret['sample'] = _sample_to_dict(seq.sample)
    ret['read_start'] = re.compile('[N\-]*').match(
        seq.sequence).span()[1] or 0

    ret['v_extent'] = ret['v_length'] + ret['num_gaps'] + ret['pad_length']
    if seq.mutations_from_clone is not None:
        ret['mutations'] = json.loads(seq.mutations_from_clone)
    else:
        ret['mutations'] = []

    ret['clone'] = _clone_to_dict(seq.clone) if seq.clone is not None else None
    ret['collapse_info'] = funcs.trace_seq_collapses(session, seq)

    return ret


def get_all_sequences(session, filters, order_field, order_dir, paging=None):
    """Gets a list of all clones"""
    def get_field(key):
        tbls = [Sequence, Subject, Clone]
        for t in tbls:
            if hasattr(t, key):
                return getattr(t, key)

    res = []
    query = session.query(Sequence).join(Sample).outerjoin(Clone)

    copy_number_field = 'copy_number'
    if filters is not None:
        if 'collapsed' in filters:
            if filters['collapsed'] == 'all':
                copy_number_field = Sequence.copy_number
            else:
                copy_number_field = getattr(
                    Sequence, 'copy_number_in_{}'.format(filters['collapsed'])
                )
        for key, value in filters.iteritems():
            if value in [None, True, False]:
                continue
            value = str(value).strip()
            if len(value) > 0 and value is not None:
                if key == 'sample_id':
                    query = query.filter(Sequence.sample_id == int(value))
                elif key == 'in_frame':
                    query = query.filter(Sequence.in_frame == int(value))
                elif key == 'min_copy_number':
                    query = query.filter(copy_number_field >= int(value))
                elif key == 'max_copy_number':
                    query = query.filter(copy_number_field <= int(value))
                elif key == 'collapsed':
                    query = query.filter(copy_number_field > 0)
                else:
                    query = query.filter(get_field(key).like(
                        value.replace('*', '%')))

    if (filters is None or
            'show_partials' not in filters or
            not filters['show_partials']):
        query = query.filter(Sequence.partial == 0)
    if (filters is None or 'show_indel' not in filters or
            not filters['show_indel']):
        query = query.filter(Sequence.probable_indel_or_misalign == 0)

    if order_field is not None:
        if order_field == 'copy_number':
            order_field = copy_number_field
        else:
            order_field = getattr(Sequence, order_field)

        if order_dir == 'asc':
            query = query.order_by(order_field)
        else:
            query = query.order_by(desc(order_field))

    for row in _page_query(query, paging):
        fields = _fields_to_dict(
            ['seq_id', 'paired', 'v_gene', 'j_gene', 'v_match', 'j_match',
             'v_length', 'j_length', 'cdr3_num_nts', 'cdr3_aa', 'in_frame',
             'functional', 'stop', 'partial',
             'probable_indel_or_misalign', 'copy_number',
             'copy_number_in_sample', 'copy_number_in_subject'],
            row)

        fields['sample'] = _sample_to_dict(row.sample)
        res.append(fields)

    return res
