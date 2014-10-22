import json
import math
from sqlalchemy import desc
from sqlalchemy.sql import func

import sldb.util.lookups as lookups
from sldb.common.models import *
from sldb.common.mutations import MutationType, Mutations


def get_all_clones(session, paging=None):
    res = []
    clone_q = session.query(Clone).order_by(Clone.v_gene, Clone.j_gene,
                                            Clone.cdr3)

    if paging is not None:
        page, per_page = paging
        clone_q = clone_q.offset((page - 1) * per_page).limit(per_page)

    for c in clone_q:
        stats_comb = []
        for stat in session.query(CloneFrequency)\
                           .filter(CloneFrequency.clone_id == c.id)\
                           .filter(CloneFrequency.filter_type == 'clones_all')\
                           .order_by(desc(CloneFrequency.copy_number)):
                stats_comb.append({
                    'sample': {
                        'id': stat.sample.id,
                        'name': stat.sample.name
                    },
                    'size': stat.size,
                    'copy_number': stat.copy_number
                })

        res.append({
            'id': c.id,
            'v_gene': c.v_gene,
            'j_gene': c.j_gene,
            'cdr3': c.cdr3,
            'stats': stats_comb
        })

    return res


def compare_clones(session, uids):
    clones = {}
    for clone_id, sample_id in uids:
        clone = session.query(Clone).filter(Clone.id == clone_id).first()
        mutations = Mutations(clone.germline, clone.cdr3_num_nts)
        if clone_id not in clones:
            clones[clone_id] = {
                'clone': {
                    'id': clone.id,
                    'v_gene': clone.v_gene,
                    'j_gene': clone.j_gene,
                    'cdr3': clone.cdr3,
                    'cdr3_num_nts': clone.cdr3_num_nts,
                    'germline': clone.germline
                },
                'mutation_stats': {},
                'seqs': []
            }

        for s in session.query(Sequence)\
                        .filter(Sequence.sample_id == sample_id)\
                        .filter(Sequence.clone_id == clone_id):
            clones[clone_id]['seqs'].append({
                'seq_id': s.seq_id,
                'sample': {
                    'id': s.sample.id,
                    'name': s.sample.name
                },
                'junction_nt': s.junction_nt,
                'sequence': s.sequence,
                'mutations': mutations.add_sequence(s.sequence),
            })

        region_stats, pos_stats = mutations.get_aggregate()
        clones[clone_id]['mutation_stats']['regions'] = region_stats
        clones[clone_id]['mutation_stats']['positions'] = pos_stats

    return clones


def get_clone_overlap(session, filter_type, samples, paging=None):
    res = []
    q = session.query(CloneFrequency,
                      func.sum(CloneFrequency.copy_number).label('cn'),
                      func.group_concat(CloneFrequency.sample_id)
                      .label('samples'))\
        .filter(CloneFrequency.sample_id.in_(samples))\
        .filter(CloneFrequency.filter_type == filter_type)\
        .group_by(CloneFrequency.clone_id)\
        .order_by(desc(func.sum(CloneFrequency.copy_number)))

    if paging is not None:
        page, per_page = paging
        q = q.offset((page - 1) * per_page).limit(per_page)
        num_pages = math.ceil(len(q.all()) / per_page)

    for r in q:
        cn = int(r.cn)
        samples = ','.join(map(str, sorted(map(int, r.samples.split(',')))))
        r = r.CloneFrequency
        res.append({
            'samples': samples,
            'copy_number': cn,
            'clone': {
                'id': r.clone.id,
                'v_gene': r.clone.v_gene,
                'j_gene': r.clone.j_gene,
                'cdr3': r.clone.cdr3
            },
        })

    if paging:
        return res, num_pages
    return res


def get_v_usage(session, filter_type, samples):
    data = {}
    headers = []
    for s in session.query(SampleStats)\
            .filter(SampleStats.filter_type == filter_type)\
            .filter(SampleStats.sample_id.in_(samples)):
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
