import datetime
import json
import math
import re
import numpy as np
from sqlalchemy.sql.expression import false, true
from sqlalchemy import desc, distinct, inspect
from sqlalchemy.sql import func
from sqlalchemy.sql.expression import false
from sqlalchemy.ext.declarative import DeclarativeMeta

import sldb.util.lookups as lookups
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
    d = _fields_to_dict(['id', 'cdr3_nt'], clone)
    d['group'] = {
        'id': clone.group.id,
        'v_gene': clone.group.v_gene,
        'j_gene': clone.group.j_gene,
        'cdr3_aa': clone.group.cdr3_aa,
        'cdr3_num_nts': clone.group.cdr3_num_nts,
        'subject': _subject_to_dict(clone.group.subject),
    }
    d['germline'] = clone.group.germline[0:VGene.CDR3_OFFSET] + \
        clone.cdr3_nt + clone.group.germline[VGene.CDR3_OFFSET +
                                             clone.group.cdr3_num_nts:]
    return d


def get_all_studies(session):
    result = {}
    for sample in session.query(Sample).order_by(Sample.date):
        if session.query(Sequence).filter(
                Sequence.sample == sample).count() > 0:
            status = 'reads'
        elif session.query(NoResult).filter(
                NoResult.sample == sample).count() > 0:
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
        for stat in session.query(CloneStats).filter(
                CloneStats.clone_id == c.id,
                CloneStats.sample_id != 0
                ).order_by(desc(CloneStats.unique_cnt)):
            stats_comb.append({
                'sample': {
                    'id': stat.sample.id,
                    'name': stat.sample.name
                },
                'unique_sequences': int(stat.unique_cnt),
                'total_sequences': int(stat.total_cnt)
            })
        clone_dict = _clone_to_dict(c)
        totals = session.query(
            CloneStats.total_cnt, CloneStats.unique_cnt
        ).filter(
            CloneStats.clone_id == c.id,
            CloneStats.sample_id == 0
        ).first()
        clone_dict['unique_sequences'] = totals.unique_cnt
        clone_dict['total_sequences'] = totals.total_cnt
        clone_dict['stats'] = stats_comb
        res.append(clone_dict)

    return res


def compare_clones(session, uids):
    """Compares sequences within clones by determining their mutations"""

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
    clones = {}
    for clone_id, sample_ids in uids.iteritems():
        full_clone = None in sample_ids
        clone = session.query(Clone).filter(Clone.id == clone_id).first()

        q = session.query(distinct(Sequence.sequence_replaced)).filter(
            Sequence.clone_id == clone_id).count()
        q = session.query(
            Sequence,
            func.sum(Sequence.copy_number).label('copy_number'))\
            .filter(Sequence.clone_id == clone_id)
        if not full_clone:
            q = q.filter(Sequence.sample_id.in_(sample_ids))
        q = q.order_by(desc('copy_number')).group_by(
            Sequence.sequence_replaced)

        clones[clone_id] = {
            'clone': _clone_to_dict(clone),
            'seqs': []
        }

        start_ptrn = re.compile('[N\-]*')
        for seqr in q:
            seq = seqr.Sequence
            read_start = start_ptrn.match(seq.sequence).span()[1] or 0
            clones[clone_id]['seqs'].append({
                'seq_id': seq.seq_id,
                'sample': {
                    'id': seq.sample.id,
                    'name': seq.sample.name,
                },
                'junction_nt': seq.junction_nt,
                'sequence': seq.sequence_replaced,
                'read_start': read_start,
                'copy_number': int(seqr.copy_number),
                'mutations': json.loads(seq.mutations_from_clone),
                'v_extent': seq.v_length + seq.num_gaps + seq.pad_length,
                'j_length': seq.j_length,
            })

        if full_clone:
            full_clone_stats = session.query(
                CloneStats.mutations,
                CloneStats.unique_cnt
            ).filter(
                CloneStats.clone_id == clone_id,
                CloneStats.sample_id == 0
            ).first()
            total_seqs = full_clone_stats.unique_cnt
            all_mutations = json.loads(full_clone_stats.mutations)
        else:
            all_mutations = CloneMutations(session, clone).calculate(
                limit_samples=sample_ids, only_clone=True).get_all()
            total_seqs = session.query(func.count(Sequence.seq_id)).filter(
                Sequence.clone_id == clone_id,
                Sequence.sample_id.in_(sample_ids)).scalar()

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
        clones[clone_id]['mutation_stats'] = mut_dict

    return clones


def get_clone_tree(session, clone_id):
    return session.query(Clone.tree).filter(Clone.id == clone_id).first()


def get_clone_overlap(session, filter_type, ctype, limit,
                      paging=None):
    """Gets a list of clones and the samples in `samples` which they appear"""
    fltr = _clone_filters[filter_type]
    res = []
    q = fltr(session.query(
        CloneStats,
        func.sum(CloneStats.unique_cnt).label('unique'),
        func.sum(CloneStats.total_cnt).label('total'),
    ).join(Clone))

    if ctype == 'samples':
        q = q.filter(CloneStats.sample_id.in_(limit))
    elif ctype == 'subject':
        q = q.filter(CloneStats.sample.has(subject_id=limit))

    q = q.group_by(CloneStats.clone_id).order_by(desc('total'))

    if paging is not None:
        page, per_page = paging
        q = q.offset((page - 1) * per_page).limit(per_page)

    for clone in q:
        selected_samples = []
        other_samples = []
        for stat in session.query(CloneStats).filter(
                CloneStats.clone_id == clone.CloneStats.clone_id).order_by(
                desc(CloneStats.total_cnt)):
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
            'unique_sequences': int(clone.unique),
            'total_sequences': int(clone.total),
            'clone': _clone_to_dict(clone.CloneStats.clone),
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
        name_func = lambda s: s.split('*')[0].split('-', 1)[0].replace(
                        'IGHV', '')
    else:
        name_func = lambda s: s.split('*')[0].replace(
                        'IGHV', '')
    data = {}
    totals = {}
    headers = []
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
            totals[group_key] = 0
        for v in dist:
            totals[group_key] += v[1]

        for v in dist:
            name, occ = v
            name = '|'.join(sorted(set(map(name_func, name.split('|')))))

            if name not in data[group_key]:
                data[group_key][name] = 0
            data[group_key][name] += round(100 * occ / float(totals[group_key]), 2)

            if data[group_key][name] >= 1.0 and name not in headers:
                headers.append(name)

    for s, vs in data.iteritems():
        for header in headers:
            if header not in vs:
                vs[header] = 0

    return data, sorted(headers), totals


def get_v_usage_grouped(session, samples, filter_type, include_outliers,
                include_partials, grouping):
    """Gets the V-Gene usage percentages for samples"""
    data = {}
    groups = {}
    headers = []
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
        if group_key not in groups:
            groups[group_key] = []
        groups[group_key].append(s.sample.name)

        data[s.sample.name] = {}
        total = 0
        for v in dist:
            total += v[1]

        for v in dist:
            name, occ = v
            name = '|'.join(
                sorted(set(map(
                    lambda s: s.split('*')[0].replace(
                        'IGHV', ''), name.split('|')))))

            if name not in data[s.sample.name]:
                data[s.sample.name][name] = 0
            data[s.sample.name][name] += round(100 * occ / float(total), 2)

            if data[s.sample.name][name] >= 1.0 and name not in headers:
                headers.append(name)

    for s, vs in data.iteritems():
        for header in headers:
            if header not in vs:
                vs[header] = 'none'

    keys = []
    lookup = {}
    for i, (group, members) in enumerate(groups.iteritems()):
        for member in sorted(members):
            keys.append(member)
            lookup[member] = i
    return data, sorted(headers), keys, lookup


def get_all_subjects(session, paging):
    q = session.query(Subject)

    if paging is not None:
        page, per_page = paging
        q = q.offset((page - 1) * per_page).limit(per_page)

    subjects = []
    for subject in q:
        seqs = session.query(func.sum(SampleStats.sequence_cnt))\
            .filter(
                SampleStats.sample.has(subject=subject),
                SampleStats.filter_type == 'all',
                SampleStats.outliers == true(),
                SampleStats.full_reads == false()).scalar()
        if seqs is None:
            continue

        subjects.append({
            'id': subject.id,
            'identifier': subject.identifier,
            'study': {
                'id': subject.study.id,
                'name': subject.study.name
            },
            'total_samples': session.query(func.count(Sample.id)).filter(
                Sample.subject == subject).scalar(),
            'unique_seqs': int(seqs),
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
        'v_match_dist', 'v_length_dist', 'j_match_dist',
        'j_length_dist', 'v_gene_dist', 'j_gene_dist',
        'cdr3_length_dist', 'copy_number_dist'
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
        'seq_id', 'alignment', 'v_gene', 'j_gene',
        'junction_nt', 'junction_aa', 'germline', 'v_match', 'j_match',
        'v_length', 'j_length', 'in_frame', 'functional', 'stop',
        'copy_number', 'sequence', 'pre_cdr3_length', 'pre_cdr3_match',
        'post_cdr3_length', 'post_cdr3_match', 'pad_length', 'num_gaps',
        'probable_indel_or_misalign'], seq)
    ret['sample'] = _sample_to_dict(seq.sample)

    ret['v_extent'] = ret['v_length'] + ret['num_gaps'] + ret['pad_length']

    if seq.clone is None:
        ret['clone'] = None
    else:
        ret['clone'] = _clone_to_dict(seq.clone)

    ret['possible_duplicates'] = []
    ret['total_copy_number'] = ret['copy_number']
    for dup in session.query(Sequence).filter(
            Sequence.sequence_replaced == seq.sequence_replaced,
            Sequence.sample.has(subject_id=seq.sample.subject_id),
            Sequence.seq_id != seq.seq_id)\
            .order_by(Sequence.sample_id):
        ret['possible_duplicates'].append({
            'seq_id': dup.seq_id,
            'sample': {
                'id': dup.sample.id,
                'name': dup.sample.name,
            },
            'alignment': dup.alignment,
            'copy_number': dup.copy_number,
        })
        ret['total_copy_number'] += dup.copy_number

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

    if filters is not None:
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
                    query = query.filter(Sequence.copy_number >= int(value))
                elif key == 'max_copy_number':
                    query = query.filter(Sequence.copy_number <= int(value))
                else:
                    query = query.filter(get_field(key).like(
                        value.replace('*', '%')))

    if (filters is None
            or 'show_partials' not in filters
            or not filters['show_partials']):
        query = query.filter(Sequence.alignment == 'R1+R2')
    if (filters is None or 'show_indel' not in filters
            or not filters['show_indel']):
        query = query.filter(Sequence.probable_indel_or_misalign == 0)

    if paging is not None:
        page, per_page = paging
        query = query.offset((page - 1) * per_page).limit(per_page)

    for row in query:
        fields = _fields_to_dict(
            ['seq_id', 'alignment', 'v_gene', 'j_gene', 'v_match', 'j_match',
             'v_length', 'j_length', 'junction_num_nts', 'junction_aa',
             'in_frame', 'functional', 'stop', 'probable_indel_or_misalign'],
            row)

        fields['copy_number'] = int(session.query(
            func.sum(Sequence.copy_number)).filter(
                Sequence.unique_id == row.unique_id).scalar())

        fields['sample'] = _sample_to_dict(row.sample)
        res.append(fields)

    return res
