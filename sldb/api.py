import argparse
import json
import math

from flask import Flask, Response, request, jsonify
import flask.ext.sqlalchemy
import flask.ext.restless

from sqlalchemy import create_engine, desc
from sqlalchemy.orm import sessionmaker, scoped_session
from sqlalchemy.sql import func

from sldb.models import *
import sldb.util.lookups

app = flask.Flask(__name__)


def _add_cors_header(response):
    response.headers['Access-Control-Allow-Origin'] = '*'
    response.headers['Access-Control-Allow-Methods'] = 'GET'
    response.headers['Access-Control-Allow-Headers'] = (
        'Origin, X-Requested-With, Content-Type, Accept')
    response.headers['Access-Control-Allow-Credentials'] = 'true'

    return response


def _get_paging():
    page = request.args.get('page') or 1
    per_page = request.args.get('per_page') or 10
    page = int(page)
    per_page = int(per_page)
    return page, per_page


class MutationType(object):
    MUT_UNK = ('?', 'unknown')
    MUT_SYN = ('S', 'synonymous')
    MUT_CONS = ('C', 'conservative')
    MUT_UNCONS = ('U', 'nonconservative')
    MUT_NONE = (' ', 'none')

    @classmethod
    def get_symbol(cls, mtype):
        return mtype[0]

    @classmethod
    def get_readable(cls, mtype):
        return mtype[1]

    @classmethod
    def get_types(cls):
        return [getattr(MutationType, attr) for attr in filter(lambda a:
                a.startswith('MUT_'), dir(MutationType))]


class Mutations(object):
    def __init__(self, germline, cdr3_num_nts):
        self.germline = germline
        self.region_stats = {}
        for region in ['all', 'CDR1', 'CDR2', 'CDR3', 'FR1', 'FR2', 'FR3']:
            self.region_stats[region] = self._create_count_record()
        self.pos_stats = {}
        self.cdr3_num_nts = cdr3_num_nts

    def _get_region(self, index):
        if index <= 77:
            return 'FR1'
        elif index <= 113:
            return 'CDR1'
        elif index <= 164:
            return 'FR2'
        elif index <= 194:
            return 'CDR2'
        elif index <= 308:
            return 'FR3'
        elif index <= 308 + self.cdr3_num_nts:
            return 'CDR3'
        return 'FR4'

    def _create_count_record(self, int_count=False):
        rec = {}
        for m in (MutationType.MUT_SYN, MutationType.MUT_CONS,
                  MutationType.MUT_UNCONS):
            if int_count:
                rec[MutationType.get_readable(m)] = 0
            else:
                rec[MutationType.get_readable(m)] = []
        return rec

    def _add_region_stat(self, i, seq):
        mtype = self._get_mut_type(seq, i)
        if mtype not in (MutationType.MUT_NONE, MutationType.MUT_UNK):
            mtype = MutationType.get_readable(mtype)
            region = self._get_region(i)
            if region not in self.region_stats:
                self.region_stats[region] = self._create_count_record()

            mutation = (i, self.germline[i], seq[i])
            self.region_stats[region][mtype].append(mutation)
            self.region_stats['all'][mtype].append(mutation)

    def _add_pos_stat(self, i, mtype, seq):
        mtype = self._get_mut_type(seq, i)
        if mtype not in (MutationType.MUT_NONE, MutationType.MUT_UNK):
            if i not in self.pos_stats:
                self.pos_stats[i] = self._create_count_record(int_count=True)
            self.pos_stats[i][MutationType.get_readable(mtype)] += 1

    def _get_aa_at(self, seq, i):
        aa_off = i - i % 3
        return aa_from_codon(seq[aa_off:aa_off + 3])

    def _get_mut_type(self, seq, i):
        if self.germline[i] != seq[i]:
            grm_aa = self._get_aa_at(self.germline, i)
            seq_aa = self._get_aa_at(seq, i)

            if grm_aa is None or seq_aa is None:
                return MutationType.MUT_UNK
            elif grm_aa != seq_aa:
                if lookups.are_conserved_aas(grm_aa, seq_aa):
                    return MutationType.MUT_CONS
                return MutationType.MUT_UNCONS
            else:
                return MutationType.MUT_SYN
        else:
            return MutationType.MUT_NONE

    def add_sequence(self, seq):
        mut_str = ''
        for i in range(0, len(self.germline)):
            mut = self._get_mut_type(seq, i)
            mut_str += MutationType.get_symbol(mut)
            self._add_region_stat(i, seq)
            self._add_pos_stat(i, mut, seq)
        return mut_str

    def get_aggregate(self):
        final_region_stats = {}
        for r, regions in self.region_stats.iteritems():
            final_region_stats[r] = {
                'counts': {
                    'unique': self._create_count_record(True),
                    'total': self._create_count_record(True)
                },
                'mutations': {}
            }

            for mtype, stats in regions.iteritems():
                final_region_stats[r]['mutations'][mtype] = []
                for mutation in stats:
                    count = len(filter(lambda e: e == mutation, stats))
                    mutation_type_cnts = \
                        final_region_stats[r]['mutations'][mtype]
                    if (count, mutation) not in mutation_type_cnts:
                        mutation_type_cnts.append((count, mutation))

        for r, region in final_region_stats.iteritems():
            for mtype, stats in region['mutations'].iteritems():
                region['counts']['total'][mtype] = reduce(
                    lambda a, b: a + b[0],
                    region['mutations'][mtype], 0)
                region['counts']['unique'][mtype] = \
                    len(region['mutations'][mtype])

            region['counts']['total']['sum'] = \
                sum(region['counts']['total'].values())
            region['counts']['unique']['sum'] = \
                sum(region['counts']['unique'].values())

        return final_region_stats, self.pos_stats


@app.route('/api/clones/', methods=['GET'])
def clones():
    session = scoped_session(session_factory)()
    page, per_page = _get_paging()

    res = []
    for c in session.query(Clone).order_by(Clone.v_gene, Clone.j_gene,
                                           Clone.cdr3)\
            .offset((page - 1) * per_page)\
            .limit(per_page):

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

    session.close()
    return jsonify(objects=res)


@app.route('/api/clone_compare/<uids>', methods=['GET'])
def clone_compare(uids):
    session = scoped_session(session_factory)()
    clones = {}
    for u in uids.split(','):
        clone_id, sample_id = map(int, u.split('_'))

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

    session.close()
    return jsonify(clones=clones)


def _clone_overlap_data(filter_type, samples, download):
    session = scoped_session(session_factory)()
    if not download:
        page, per_page = _get_paging()

    samples = map(int, samples.split(','))
    res = []
    q = session.query(CloneFrequency,
                      func.sum(CloneFrequency.copy_number).label('cn'),
                      func.group_concat(CloneFrequency.sample_id)
                      .label('samples'))\
        .filter(CloneFrequency.sample_id.in_(samples))\
        .filter(CloneFrequency.filter_type == filter_type)\
        .group_by(CloneFrequency.clone_id)\
        .order_by(desc(func.sum(CloneFrequency.copy_number)))
    num_pages = math.ceil(len(q.all()) / per_page)
    if not download:
        q = q.offset((page - 1) * per_page).limit(per_page)
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

    session.close()

    return res, num_pages


@app.route('/api/clone_overlap/<filter_type>/<samples>', methods=['GET'])
def clone_overlap(filter_type, samples):
    items, num_pages = _clone_overlap_data(filter_type, samples, False)
    return jsonify(items=items, num_pages=num_pages)


@app.route('/api/data/clone_overlap/<filter_type>/<samples>', methods=['GET'])
def download_clone_overlap(filter_type, samples):
    data, _ = _clone_overlap_data(filter_type, samples, True)

    def _gen(data):
        yield ','.join(['samples', 'copy_number', 'v_gene', 'j_gene',
                       'cdr3']) + '\n'
        for c in data:
            yield ','.join(map(str, [c['samples'].replace(',', ' '),
                           c['copy_number'],
                           c['clone']['v_gene'],
                           c['clone']['j_gene'],
                           c['clone']['cdr3']])) + '\n'

    return Response(_gen(data), headers={
        'Content-Disposition':
        'attachment;filename={}_{}.csv'.format(
            filter_type,
            samples.replace(',', '-'))})


def _v_usage_data(filter_type, samples):
    session = scoped_session(session_factory)()
    data = {}
    samples = map(int, samples.split(','))
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

    session.close()
    return data, headers


@app.route('/api/v_usage/<filter_type>/<samples>', methods=['GET'])
def v_usage(filter_type, samples):
    data, headers = _v_usage_data(filter_type, samples)
    x_categories = headers
    y_categories = data.keys()

    array = []
    for j, y in enumerate(y_categories):
        s = 0
        for i, x in enumerate(x_categories):
            usage_for_y = data[y]
            if x in usage_for_y:
                array.append([i, j, usage_for_y[x]])
                s += usage_for_y[x]
            else:
                array.append([i, j, 0])

    return jsonify(x_categories=x_categories,
                   y_categories=y_categories,
                   data=array)


@app.route('/api/data/v_usage/<filter_type>/<samples>', methods=['GET'])
def download_v_usage(filter_type, samples):
    data, headers = _v_usage_data(filter_type, samples)
    ret = 'sample,' + ','.join(headers) + '\n'
    for sample, dist in data.iteritems():
        row = [sample]
        for gene in headers:
            if gene in dist:
                row.append(dist[gene])
            else:
                row.append(0)
        ret += ','.join(map(str, row)) + '\n'

    return Response(ret, headers={
        'Content-Disposition':
        'attachment;filename=v_usage-{}_{}.csv'.format(
            filter_type,
            ','.join(map(str, samples)).replace(',', '-'))})


def _init_db(user, pw, db):
    engine = create_engine(('mysql://{}:{}@localhost/'
                            '{}?charset=utf8&use_unicode=0').format(
                                user, pw, db))

    Base.metadata.create_all(engine)
    Base.metadata.bind = engine
    return sessionmaker(bind=engine)


def run_api():
    parser = argparse.ArgumentParser(
        description='Provides a restless interface to the master table '
                    'database')
    parser.add_argument('db', help='mySQL database')
    parser.add_argument('user', help='mySQL user')
    parser.add_argument('pw', help='mySQL password')
    parser.add_argument('-p', default=5000, type=int, help='API offer port')
    args = parser.parse_args()

    app.config['SQLALCHEMY_DATABASE_URI'] = (
        'mysql://{}:{}@localhost/'
        '{}?charset=utf8&use_unicode=0').format(args.user, args.pw, args.db)
    db = flask.ext.sqlalchemy.SQLAlchemy(app)

    manager = flask.ext.restless.APIManager(app, flask_sqlalchemy_db=db)

    manager.create_api(Study, methods=['GET'])
    manager.create_api(Sample, methods=['GET'], include_columns=[
                       'id', 'name', 'info', 'study'])
    manager.create_api(Sequence, methods=['GET'])
    manager.create_api(SampleStats, methods=['GET'], collection_name='stats',
                       max_results_per_page=10000,
                       results_per_page=10000)
    manager.create_api(CloneFrequency, methods=['GET'],
                       collection_name='clone_freqs',
                       exclude_columns=['sample'])
    session_factory = _init_db(args.user, args.pw, args.db)

    app.after_request(_add_cors_header)
    app.run(host='0.0.0.0', port=args.p, debug=True, threaded=True)
