import argparse
import json
import math

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, scoped_session

import bottle
from bottle import route, response, request, install, run

import sldb.api.queries as queries
from sldb.common.models import *


class EnableCors(object):
    name = 'enable_cors'
    api = 2

    def apply(self, fn, context):
        def _enable_cors(*args, **kwargs):
            response.headers['Access-Control-Allow-Origin'] = '*'
            response.headers['Access-Control-Allow-Methods'] = ('GET, POST')
            response.headers['Access-Control-Allow-Headers'] = (
                'Origin, '
                'Accept, Content-Type, '
                'X-Requested-With, X-CSRF-Token')

            if bottle.request.method != 'OPTIONS':
                return fn(*args, **kwargs)

        return _enable_cors


def _get_arg(key, is_json=True):
    if key not in request.query or len(request.query[key].strip()) == 0:
        return None
    req = request.query[key].strip()
    return json.loads(req) if is_json else req


def _get_paging():
    """Handles paging based on a request's query string"""
    page = _get_arg('page', False) or 1
    per_page = _get_arg('per_page', False) or 10
    page = int(page)
    per_page = int(per_page)
    return page, per_page


def _split(ids, delim=','):
    """Helper function to split a string into an integer array"""
    return map(int, ids.split(delim))


@route('/api/sequences/')
def sequences():
    """Gets a list of all sequences"""
    session = scoped_session(session_factory)()
    sequences = queries.get_all_sequences(
        session,
        _get_arg('filter'),
        _get_arg('order_field', False) or 'seq_id',
        _get_arg('order_dir', False) or 'desc',
        _get_paging())
    session.close()
    return json.dumps({'sequences': sequences})


@route('/api/sequence/<sample_id>/<seq_id>')
def sequence(sample_id, seq_id):
    session = scoped_session(session_factory)()
    seq = queries.get_sequence(session, int(sample_id), seq_id)
    session.close()
    return json.dumps({'sequence': seq})


@route('/api/studies')
def studies():
    session = scoped_session(session_factory)()
    studies = queries.get_all_studies(session)
    session.close()
    return json.dumps({'studies': studies})


@route('/api/subjects')
def subjects():
    """Gets a list of all subjects"""
    session = scoped_session(session_factory)()
    subjects = queries.get_all_subjects(session, _get_paging())
    session.close()
    return json.dumps({'subjects': subjects})


@route('/api/subject/<sid>')
def subject(sid):
    session = scoped_session(session_factory)()
    subject = queries.get_subject(session, int(sid))
    session.close()
    return json.dumps({'subject': subject})


@route('/api/clones/')
def clones():
    """Gets a list of all clones"""
    session = scoped_session(session_factory)()
    clones = queries.get_all_clones(
        session,
        _get_arg('filter'),
        _get_arg('order_field', False) or 'id',
        _get_arg('order_dir', False) or 'desc',
        _get_paging())
    session.close()
    return json.dumps({'clones': clones})


@route('/api/clone_compare/<uids>')
def clone_compare(uids):
    """Compares clones by determining their mutations"""
    session = scoped_session(session_factory)()
    clones_and_samples = {}
    for u in uids.split(','):
        if u.find('_') < 0:
            clone = int(u)
            sample = None
        else:
            clone, sample = _split(u, '_')
        if clone not in clones_and_samples:
            clones_and_samples[clone] = set([])
        clones_and_samples[clone].add(sample)

    clones = queries.compare_clones(session, clones_and_samples)
    session.close()
    return json.dumps({'clones': clones})


@route('/api/clone_overlap/<filter_type>/<samples>')
@route('/api/subject_clones/<filter_type>/<subject>')
def clone_overlap(filter_type, samples=None, subject=None):
    """Gets clonal overlap between samples"""
    session = scoped_session(session_factory)()
    if samples is not None:
        sids = _split(samples)
    ctype = 'samples' if samples is not None else 'subject'

    clones = queries.get_clone_overlap(
        session, filter_type, ctype,
        sids if samples is not None else subject, _get_paging())
    session.close()
    return json.dumps({'clones': clones})


@route('/api/stats/<samples>')
def stats(samples):
    session = scoped_session(session_factory)()
    stats = queries.get_stats(session, _split(samples), False)
    outliers = queries.get_stats(session, _split(samples), True)
    session.close()
    return json.dumps({'no_outliers': stats, 'with_outliers': outliers})


@route('/api/data/clone_overlap/<filter_type>/<samples>')
@route('/api/data/subject_clones/<filter_type>/<subject>')
def download_clone_overlap(filter_type, samples=None, subject=None):
    """Downloads a CSV of the clonal overlap between samples"""
    if samples is not None:
        ctype = 'samples'
        limit = _split(samples) if samples is not None else None
        fn = 'overlap_{}_{}.csv'.format(
            filter_type,
            samples.replace(',', '-'))
    else:
        ctype = 'subject'
        limit = int(subject)
        fn = 'subject_{}_{}.csv'.format(
            filter_type,
            subject)

    session = scoped_session(session_factory)()
    data = queries.get_clone_overlap(session, filter_type, ctype, limit)
    session.close()

    response.headers['Content-Disposition'] = \
        'attachment;filename={}'.format(fn)
    yield ','.join([
        'clone_id', 'samples', 'total_sequences', 'unique_sequences',
        'subject', 'v_gene', 'j_gene', 'cdr3_len', 'cdr3_aa', 'cdr3_nt']) + \
        '\n'

    for c in data:
        fields = [
            c['clone']['id'],
            ' '.join(c['samples']),
            c['total_sequences'],
            c['unique_sequences'],
            '{} ({})'.format(
                c['clone']['group']['subject']['identifier'],
                c['clone']['group']['subject']['study']['name']),
            c['clone']['group']['v_gene'],
            c['clone']['group']['j_gene'],
            c['clone']['group']['cdr3_aa'],
            c['clone']['cdr3_nt']]

        yield ','.join(map(str, fields)) + '\n'


@route(
    '/api/data/sequences/<file_type>/<replace_germ>/<cid>/<params>',
    methods=['GET'])
def download_sequences(file_type, replace_germ, cid, params):
    cid = int(cid)
    replace_germ = True if replace_germ == 'true' else False
    fn = '{}-{}'.format(cid, 'replaced' if replace_germ else 'original')

    session = scoped_session(session_factory)()
    sequences = []
    germline = None
    for clone_and_sample in params.split(','):
        if '_' in clone_and_sample:
            clone_id, sample = _split(clone_and_sample, '_')
        else:
            clone_id = int(clone_and_sample)
            sample = None
        if clone_id != cid:
            continue
        if sample is not None:
            fn += '_{}'.format(sample)
        else:
            fn += '_all'

        query = session.query(SequenceMapping).filter(
            SequenceMapping.clone_id == cid)
        if sample is not None:
            query = query.filter(SequenceMapping.sample_id == sample)
        for s in query:
            if germline is None:
                germline = s.clone.group.germline
            sequences.append({
                'sample_name': s.sample.name,
                'seq_id': s.seq_id,
                'sequence': s.identity_seq.sequence_replaced
                if replace_germ else s.sequence
            })
    session.close()

    def _gen(sequences):
        yield '>>germline\n{}\n\n'.format(germline)
        for seq in sequences:
            if file_type == 'fasta':
                yield '>{},{}\n{}\n\n'.format(
                    seq['sample_name'],
                    seq['seq_id'],
                    seq['sequence'])
    response.headers['Content-Disposition'] = \
        'attachment;filename={}.fasta'.format(fn)
    return _gen(sequences)


@route('/api/v_usage/<filter_type>/<samples>')
def v_usage(filter_type, samples):
    """Gets the V usage for samples in a heatmap-formatted array"""
    session = scoped_session(session_factory)()
    data, headers = queries.get_v_usage(session, filter_type, _split(samples))
    session.close()
    x_categories = headers
    y_categories = data.keys()

    array = []
    for j, y in enumerate(y_categories):
        for i, x in enumerate(x_categories):
            usage_for_y = data[y]
            if x in usage_for_y:
                array.append([i, j, usage_for_y[x]])
            else:
                array.append([i, j, 0])

    return json.dumps({
        'x_categories': x_categories,
        'y_categories': y_categories,
        'data': array
    })


@route('/api/data/v_usage/<filter_type>/<samples>')
def download_v_usage(filter_type, samples):
    """Downloads a CSV of the V usage for samples"""
    session = scoped_session(session_factory)()
    data, headers = queries.get_v_usage(session, filter_type, _split(samples))
    session.close()
    ret = 'sample,' + ','.join(headers) + '\n'
    for sample, dist in data.iteritems():
        row = [sample]
        for gene in headers:
            if gene in dist:
                row.append(dist[gene])
            else:
                row.append(0)
        ret += ','.join(map(str, row)) + '\n'

    response.headers['Content-Disposition'] = \
        'attachment;filename=v_usage-{}_{}.csv'.format(
            filter_type,
            samples.replace(',', '-'))
    return ret


def run_rest_service(session_maker, args):
    """Runs the rest service based on command line arguments"""
    global session_factory
    session_factory = session_maker
    bottle.install(EnableCors())
    if args.debug:
        bottle.run(host='0.0.0.0', port=args.port, debug=True)
    else:
        bottle.run(host='0.0.0.0', port=args.port, server='gevent')
