from collections import Counter
import json
import math
import re

from sqlalchemy import and_, desc, distinct
from sqlalchemy.orm.strategy_options import Load
from sqlalchemy.sql import func
from sqlalchemy.sql.expression import false, true

from immunedb.common.models import (Clone, CloneStats, Sample, SampleStats,
                                    SelectionPressure, Sequence,
                                    SequenceCollapse, Subject)
from immunedb.common.mutations import threshold_mutations


_clone_filters = {
    'clones_all': lambda q: q,
    'clones_functional': lambda q: q.filter(Clone.functional == 1),
    'clones_nonfunctional': lambda q: q.filter(Clone.functional == 0)
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
    d = _fields_to_dict(['id', 'name'], sample)
    d['subject'] = _subject_to_dict(sample.subject)
    d['metadata'] = sample.metadata_dict
    return d


def _clone_to_dict(clone):
    d = _fields_to_dict(['id', 'cdr3_nt', 'v_gene', 'j_gene', 'cdr3_aa',
                         'cdr3_num_nts', 'regions', 'insertions', 'deletions',
                         'overall_unique_cnt', 'overall_total_cnt',
                         'overall_instance_cnt',
                         'overall_unique_cnt_with_subclones',
                         'overall_instance_cnt_with_subclones',
                         'overall_total_cnt_with_subclones'], clone)
    d['subject'] = _subject_to_dict(clone.subject)
    d['germline'] = clone.consensus_germline

    return d


def _page_query(q, paging):
    if paging is None:
        return q
    page, per_page = paging
    return q.offset(max(0, (page - 1) * per_page)).limit(per_page)


def get_samples(session, sample_ids=None):
    query = session.query(
        Sample, SampleStats
    ).outerjoin(
        SampleStats,
        and_(
            SampleStats.sample_id == Sample.id,
            SampleStats.outliers == true(),
            SampleStats.full_reads == false(),
            SampleStats.filter_type == 'all'
        )
    )

    if sample_ids:
        query = query.filter(SampleStats.sample_id.in_(sample_ids))

    query = query.order_by(Sample.subject_id).options(
        Load(SampleStats).load_only(
            'sequence_cnt', 'in_frame_cnt', 'stop_cnt', 'functional_cnt',
            'no_result_cnt'
        )
    )

    clone_counts = session.query(
        SampleStats.sample_id,
        SampleStats.sequence_cnt,
        SampleStats.functional_cnt
    ).filter(
        SampleStats.outliers == true(),
        SampleStats.full_reads == false(),
        SampleStats.filter_type == 'clones_all'
    )
    clone_counts = {s.sample_id: s for s in clone_counts}

    result = []
    for sample, stats in query:
        sample_dict = _sample_to_dict(sample)
        if stats is not None:
            clones = clone_counts[sample.id]
            sample_dict['sequence_cnt'] = stats.sequence_cnt
            sample_dict['functional_cnt'] = stats.functional_cnt
            sample_dict['no_result_cnt'] = stats.no_result_cnt
            sample_dict['total_cnt'] = stats.sequence_cnt + stats.no_result_cnt
            sample_dict['clone_cnt'] = clones.sequence_cnt
            sample_dict['functional_clone_cnt'] = clones.functional_cnt
        result.append(sample_dict)

    return result


def get_clones(session, filters, order_field, order_dir, subject_limit=None,
               paging=None):
    """Gets a list of all clones"""
    res = []
    clone_q = session.query(Clone)

    if subject_limit is not None:
        clone_q = clone_q.filter(Clone.subject_id == subject_limit)

    if filters is not None:
        size_field = filters.pop('size_field', 'copies')
        if size_field == 'copies':
            size_field = Clone.overall_total_cnt
        elif size_field == 'uniques':
            size_field = Clone.overall_unique_cnt
        elif size_field == 'instances':
            size_field = Clone.overall_instance_cnt

        for key, value in filters.items():
            if value is None:
                continue
            value = str(value).strip()
            if len(value) > 0 and value is not None:
                if key == 'min_cdr3_num_nts':
                    clone_q = clone_q.filter(Clone.cdr3_num_nts >= int(value))
                elif key == 'max_cdr3_num_nts':
                    clone_q = clone_q.filter(Clone.cdr3_num_nts <= int(value))
                elif key == 'cdr3_nt':
                    clone_q = clone_q.filter(Clone.cdr3_nt.like(value))
                elif key == 'cdr3_aa':
                    clone_q = clone_q.filter(Clone.cdr3_aa.like(value))
                elif key == 'min_size':
                    clone_q = clone_q.filter(size_field >= int(value))
                elif key == 'max_size':
                    clone_q = clone_q.filter(size_field <= int(value))
                elif key == 'id':
                    clone_q = clone_q.filter(Clone.id == int(value))
                elif key == 'subject_id':
                    clone_q = clone_q.filter(Clone.subject_id == int(value))
                else:
                    clone_q = clone_q.filter(
                        getattr(Clone, key).like(value))

    if order_field:
        order_field = getattr(Clone, order_field)
        if order_dir == 'asc':
            clone_q = clone_q.order_by(order_field)
        else:
            clone_q = clone_q.order_by(desc(order_field))

    for c in _page_query(clone_q, paging):
        stats_comb = []

        query = session.query(
            CloneStats.unique_cnt, CloneStats.total_cnt, Sample
        ).join(Sample).filter(
            CloneStats.clone_id == c.id,
            ~CloneStats.sample_id.is_(None)
        ).order_by(
            desc(CloneStats.unique_cnt), desc(CloneStats.total_cnt)
        )
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
        clone_dict['stats'] = stats_comb
        res.append(clone_dict)

    return res


def get_clone(session, clone_id):
    result = {}
    clone = session.query(Clone).filter(Clone.id == clone_id).first()

    stats = session.query(
        CloneStats
    ).filter(
        CloneStats.clone_id == clone.id,
    ).order_by(desc(CloneStats.unique_cnt))

    result = {
        'clone': _clone_to_dict(clone),
        'parent': _clone_to_dict(clone.parent) if clone.parent else None,
        'children': [_clone_to_dict(c) for c in clone.children],
        'samples': {
            'single': []
        }
    }

    for stat in stats:
        if stat.sample_id is None:
            result['samples']['all'] = {
                'unique': stat.unique_cnt,
                'total': stat.total_cnt
            }
        else:
            sample = _sample_to_dict(stat.sample)
            sample['unique'] = stat.unique_cnt
            sample['total'] = stat.total_cnt
            result['samples']['single'].append(sample)

    return result


def get_clone_mutations(session, clone_id, threshold_type, threshold_val,
                        sample_id=None):
    clone_stats = session.query(
        CloneStats.mutations,
        CloneStats.unique_cnt
    ).filter(
        CloneStats.clone_id == clone_id,
        CloneStats.sample_id == sample_id
    ).first()
    all_mutations = json.loads(clone_stats.mutations)
    total_seqs = clone_stats.unique_cnt

    result = {
        'positions': all_mutations['positions'],
        'regions': {},
        'total_seqs': total_seqs
    }
    if threshold_type == 'sequences':
        seq_min = threshold_val
    else:
        seq_min = int(math.ceil(threshold_val / 100.0 * total_seqs))
    result['regions'] = threshold_mutations(all_mutations, seq_min)

    return result


def get_clone_sequences(session, clone_id, get_collapse, paging):
    query = session.query(
        Sequence, SequenceCollapse.copy_number_in_subject,
        SequenceCollapse.instances_in_subject
    ).join(SequenceCollapse).filter(
        Sequence.clone_id == clone_id,
        SequenceCollapse.copy_number_in_subject > 0
    ).order_by(
        desc(SequenceCollapse.copy_number_in_subject)
    )
    query = _page_query(query, paging)

    sequences = {}
    start_ptrn = re.compile(r'[N-]*')
    for seq, copy_number_in_subject, instances_in_subject in query:
        sequences[(seq.sample.id, seq.seq_id)] = {
            'seq_id': seq.seq_id,
            'sample': {
                'id': seq.sample.id,
                'name': seq.sample.name,
            },
            'cdr3_nt': seq.cdr3_nt,
            'sequence': seq.clone_sequence,
            'read_start': start_ptrn.match(seq.sequence).span()[1] or 0,
            'copy_number_in_subject': int(copy_number_in_subject),
            'instances_in_subject': int(instances_in_subject),
            'mutations': json.loads(seq.mutations_from_clone),
            'v_extent': seq.get_v_extent(in_clone=True),
            'j_length': seq.j_length,
            'collapse_to': []
        }

    if get_collapse:
        q = session.query(
            Sequence, SequenceCollapse
        ).outerjoin(SequenceCollapse).filter(
            Sequence.clone_id == clone_id,
            SequenceCollapse.copy_number_in_subject == 0
        ).order_by(desc(Sequence.copy_number))

        for seq, collapse in q:
            key = (collapse.collapse_to_subject_sample_id,
                   collapse.collapse_to_subject_seq_id)
            if key in sequences:
                sequences[key]['collapse_to'].append({
                    'sample_id': seq.sample.id,
                    'sample_name': seq.sample.name,
                    'seq_id': seq.seq_id,
                    'copy_number_in_sample': seq.copy_number
                })

    return sorted(
        sequences.values(),
        key=lambda v: v['copy_number_in_subject'],
        reverse=True
    )


def get_selection_pressure(session, clone_id):
    query = session.query(
        SelectionPressure
    ).filter(
        SelectionPressure.clone_id == clone_id,
    )

    pressure = {}
    for row in query:
        if row.sample_id is None:
            name = 'All'
        else:
            name = session.query(Sample.name).filter(
                Sample.id == row.sample_id).first().name
        if name not in pressure:
            pressure[name] = {
                'pressure': {},
                'sample': {
                    'id': row.sample_id,
                    'name': name
                }
            }
        pressure[name]['pressure'][row.threshold] = row.to_dict()

    return pressure.values()


def get_clone_tree(session, clone_id):
    tree = session.query(Clone.tree).filter(Clone.id == clone_id).first().tree
    return json.loads(tree) if tree is not None else None


def get_clone_overlap(session, sample_ids, filter_type, order_by='total_cnt',
                      paging=None):
    """Gets a list of clones and the samples in `samples` which they appear"""
    stats = session.query(
        CloneStats.clone_id,
        CloneStats.functional,
        CloneStats.unique_cnt,
        CloneStats.total_cnt
    ).filter(
        CloneStats.sample_id.in_(sample_ids),
    )
    if filter_type != 'clones_all':
        stats = stats.filter(
            CloneStats.functional == (filter_type == 'clones_functional')
        )

    aggregated_stats = {}
    for clone_id, functional, unique_cnt, total_cnt in stats:
        aggregated_stats.setdefault(clone_id, Counter()).update({
            'unique_cnt': unique_cnt,
            'total_cnt': total_cnt
        })
    aggregated_stats = sorted(
        aggregated_stats.items(),
        key=lambda kv: kv[1][order_by], reverse=True
    )

    if paging:
        start, per_page = paging
        aggregated_stats = aggregated_stats[
            (start - 1) * per_page:start * per_page]

    res = []
    for clone_id, counts in aggregated_stats:
        selected_samples = []
        other_samples = []
        clone = session.query(Clone).filter(Clone.id == clone_id).one()
        query = session.query(CloneStats).filter(
            CloneStats.clone_id == clone_id,
            ~CloneStats.sample_id.is_(None)
        ).order_by(
            desc(getattr(CloneStats, order_by))
        )
        for stat in query:
            data = {
                'id': stat.sample_id,
                'name': stat.sample.name,
                'unique_sequences': stat.unique_cnt,
                'total_sequences': stat.total_cnt
            }
            if stat.sample_id in sample_ids:
                selected_samples.append(data)
            else:
                other_samples.append(data)

        res.append({
            'unique_sequences': counts['unique_cnt'],
            'total_sequences': counts['total_cnt'],
            'clone': {
                'id': clone.id,
                'v_gene': clone.v_gene,
                'j_gene': clone.j_gene,
                'cdr3_num_nts': clone.cdr3_num_nts,
                'cdr3_aa': clone.cdr3_aa,
                'subject': _subject_to_dict(clone.subject),
            },
            'selected_samples': selected_samples,
            'other_samples': other_samples,
        })

    return res


def get_clones_in_samples(session, samples):
    return [e.id for e in session.query(
                   distinct(Sequence.clone_id).label('id')).filter(
                   Sequence.sample_id.in_(samples))]


def get_clones_in_subject(session, subject_id):
    return [e.id for e in session.query(Clone).filter(
        Clone.subject_id == subject_id)]


def get_grouping(sample, grouping):
    if grouping == 'sample':
        return sample.name
    elif grouping == 'subject':
        return sample.subject.identifier
    else:
        return sample.metadata_dict.get(grouping, None)


def get_v_usage(session, samples, filter_type, include_outliers,
                include_partials, grouping, by_family):
    """Gets the V-Gene usage percentages for samples"""
    if by_family:
        def name_func(s, prefix):
            return s.split('*')[0].split('-', 1)[0].replace(prefix, '')
    else:
        def name_func(s, prefix):
            return s.split('*')[0].replace(prefix, '')
    data = {}
    totals = {}
    prefix = ''
    for s in session.query(SampleStats)\
            .filter(SampleStats.filter_type == filter_type,
                    SampleStats.outliers == include_outliers,
                    SampleStats.full_reads != include_partials,
                    SampleStats.sample_id.in_(samples)):
        dist = json.loads(s.v_gene_dist)

        group_key = get_grouping(s.sample, grouping)

        if group_key not in data:
            data[group_key] = {}

        for v in dist:
            name, occ = v
            prefix = name[:4]
            ties = set([name_func(n, prefix) for n in name.split('|')])
            name = '|'.join(sorted(ties))

            if name not in data[group_key]:
                data[group_key][name] = 0
            data[group_key][name] += occ

    headers = []
    for group_key, names in data.items():
        totals[group_key] = float(sum(names.values()))
        for name, value in names.items():
            percent = round(100 * value / totals[group_key], 2)
            names[name] = percent
            if name not in headers and percent >= 1.0:
                headers.append(name)

    return data, sorted(headers), totals, prefix


def get_all_subjects(session, paging=None):
    subjects = []
    for subject in _page_query(
            session.query(Subject).order_by(Subject.identifier), paging):
        total_seqs = session.query(
            func.sum(SampleStats.sequence_cnt)
        ).filter(
            SampleStats.sample.has(subject=subject),
            SampleStats.filter_type == 'all',
            SampleStats.outliers == true(),
            SampleStats.full_reads == false()
        ).scalar()
        unique_seqs = session.query(
            func.sum(SampleStats.sequence_cnt)
        ).filter(
            SampleStats.sample.has(subject=subject),
            SampleStats.filter_type == 'unique',
            SampleStats.outliers == true(),
            SampleStats.full_reads == false()
        ).scalar()
        info = {
            'id': subject.id,
            'identifier': subject.identifier,
            'study': {
                'id': subject.study.id,
                'name': subject.study.name
            },
        }
        if total_seqs is not None:
            info['total_samples'] = session.query(
                func.count(Sample.id)
            ).filter(
                Sample.subject == subject
            ).scalar()
            info['total_seqs'] = int(total_seqs)
            info['unique_seqs'] = int(unique_seqs)
            info['total_clones'] = session.query(
                func.count(Clone.id)
            ).filter(
                Clone.subject_id == subject.id
            ).scalar()
        subjects.append(info)

    return subjects


def get_subject(session, sid):
    s = session.query(Subject).filter(Subject.id == sid).first()
    samples = get_samples(
        session,
        [e.id for e in session.query(Sample.id).filter(
            Sample.subject_id == sid)])

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


def analyze_samples(session, samples, filter_type, include_outliers,
                    include_partials, percentages, grouping):
    counts = {}

    stats = {}
    sample_info = {}
    dist_fields = [
        'v_match_dist', 'v_length_dist', 'v_identity_dist', 'j_match_dist',
        'j_length_dist', 'v_gene_dist', 'j_gene_dist', 'cdr3_length_dist',
        'copy_number_dist', 'quality_dist', 'sp_fwr_dist', 'sp_cdr_dist'
    ]
    cnt_fields = [
        'sequence_cnt', 'in_frame_cnt', 'stop_cnt', 'functional_cnt',
        'no_result_cnt'
    ]

    group_sizes = {}
    # Iterate over all filter types in the samples
    for stat in session.query(SampleStats).filter(
            SampleStats.sample_id.in_(samples),
            SampleStats.outliers == include_outliers,
            SampleStats.full_reads != include_partials):
        # Update the number of sequences in each filter
        if stat.filter_type not in counts:
            counts[stat.filter_type] = 0
        counts[stat.filter_type] += stat.sequence_cnt

        # If the sample is not already in the sample_info dictionary, add it
        if stat.sample.id not in sample_info:
            sample_info[stat.sample.id] = _sample_to_dict(stat.sample)

        if stat.filter_type == 'all':
            # If the filter is for all sequences, add the count fields to the
            # sample dictionary.
            for field in cnt_fields:
                sample_info[stat.sample.id][field] = getattr(stat, field)
            sample_info[stat.sample.id]['total_cnt'] = (
                sample_info[stat.sample.id]['sequence_cnt'] +
                sample_info[stat.sample.id]['no_result_cnt']
            )
        elif stat.filter_type == 'clones_all':
            # If it's the clones_all filter, update the total number of clones
            # in the samples
            sample_info[stat.sample_id].update({
                'clones_cnt': stat.sequence_cnt,
                'clones_functional_cnt': stat.functional_cnt
            })
        elif stat.filter_type == 'unique':
            # If it's the unique filter, update the total number of unique
            # sequences in the sample
            sample_info[stat.sample_id].update({
                'unique_cnt': stat.functional_cnt
            })
        elif stat.filter_type == 'unique_multiple':
            # If it's the unique_multiple filter, update the total number of
            # unique sequences that occur multiple times in the sample
            sample_info[stat.sample_id].update({
                'unique_multiple_cnt': stat.functional_cnt
            })

        # If this is the selected filter, group and tally the statistics
        if stat.filter_type == filter_type:
            group_key = get_grouping(stat.sample, grouping)

            if group_key not in stats:
                stats[group_key] = {}
                group_sizes[group_key] = 0
            group_sizes[group_key] += 1

            fields = _fields_to_dict(dist_fields, stat)

            for field, values in fields.items():
                if field not in stats[group_key]:
                    stats[group_key][field] = {}

                for (x, freq) in json.loads(values):
                    if x not in stats[group_key][field]:
                        stats[group_key][field][x] = 0
                    stats[group_key][field][x] += freq

    for group, key_dict in stats.items():
        for key, vals in key_dict.items():
            if key == 'quality_dist':
                vals = {k: v / group_sizes[group] for k, v in vals.items()}
            reduced = []
            for x in sorted(vals.keys()):
                val = vals[x]
                if (key != 'quality_dist' and
                        percentages and sum(vals.values()) > 0):
                    val = round(100 * vals[x] / sum(vals.values()), 2)
                reduced.append((x, val))
            key_dict[key] = reduced

    clone_cnts = {s.functional: s.clones for s in session.query(
        CloneStats.functional,
        func.count(distinct(CloneStats.clone_id)).label('clones')
    ).filter(
        CloneStats.sample_id.in_(samples)
    ).group_by(CloneStats.functional)}

    counts['clones_nonfunctional'] = clone_cnts.get(False, 0)
    counts['clones_functional'] = clone_cnts.get(True, 0)
    counts['clones_all'] = clone_cnts.get(False, 0) + clone_cnts.get(True, 0)

    return {'samples': sample_info, 'counts': counts, 'stats': stats}


def trace_seq_collapses(session, seq):
    collapse_info = session.query(
        SequenceCollapse
    ).filter(
        SequenceCollapse.sample_id == seq.sample_id,
        SequenceCollapse.seq_ai == seq.ai,
    ).first()

    if collapse_info is None:
        return None
    return {
        'sample_id': collapse_info.collapse_to_subject_sample_id,
        'sample_name': collapse_info.collapse_to_seq.sample.name,
        'ai': collapse_info.collapse_to_subject_seq_ai,
        'seq_id': collapse_info.collapse_to_subject_seq_id,
        'copy_number':
            collapse_info.collapse_to_seq.collapse.copy_number_in_subject,
        'instances':
            collapse_info.collapse_to_seq.collapse.instances_in_subject,
    }


def get_sequence(session, sample_id, seq_id):
    seq = session.query(Sequence)\
        .filter(Sequence.sample_id == sample_id,
                Sequence.seq_id == seq_id).one()
    ret = _fields_to_dict([
        'seq_id', 'partial', 'rev_comp', 'v_gene', 'j_gene', 'cdr3_nt',
        'cdr3_aa', 'germline', 'v_match', 'j_match', 'v_length',
        'j_length', 'in_frame', 'functional', 'stop', 'copy_number',
        'sequence', 'pre_cdr3_length', 'pre_cdr3_match', 'post_cdr3_length',
        'post_cdr3_match', 'seq_start', 'num_gaps',
        'probable_indel_or_misalign', 'insertions', 'deletions',
        'quality', 'regions'
    ], seq)
    ret['sample'] = _sample_to_dict(seq.sample)
    ret['read_start'] = re.compile(r'[N-]*').match(
        seq.sequence).span()[1] or 0

    ret['v_extent'] = seq.get_v_extent(in_clone=False)
    ret['mutations'] = {}
    if seq.mutations_from_clone:
        ret['mutations'] = json.loads(seq.mutations_from_clone)

    ret['clone'] = _clone_to_dict(seq.clone) if seq.clone is not None else None
    ret['collapse_info'] = trace_seq_collapses(session, seq)

    return ret


def get_sequences(session, filters, order_field, order_dir, subject_id=None,
                  paging=None):
    """Gets a list of all clones"""
    res = []
    query = session.query(Sequence).outerjoin(SequenceCollapse)

    copy_number_field = 'copy_number'
    if filters is not None:
        if filters.pop('copy_type', 'sample') == 'sample':
            copy_number_field = Sequence.copy_number
        else:
            copy_number_field = SequenceCollapse.copy_number_in_subject

        for key, value in filters.items():
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
                    query = query.filter(getattr(Sequence, key).like(value))

    if subject_id is not None:
        query = query.filter(Sequence.subject_id == subject_id)

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
            ['seq_id', 'v_gene', 'j_gene', 'v_match', 'j_match',
             'v_length', 'j_length', 'cdr3_num_nts', 'cdr3_aa', 'in_frame',
             'functional', 'stop', 'partial', 'probable_indel_or_misalign',
             'copy_number'],
            row)
        if row.collapse is not None:
            fields.update({
                'copy_number_in_subject': row.collapse.copy_number_in_subject,
                'instances_in_subject': row.collapse.instances_in_subject
            })

        fields['sample'] = _sample_to_dict(row.sample)
        res.append(fields)

    return res
