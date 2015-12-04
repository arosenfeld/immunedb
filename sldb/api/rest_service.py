import json
from functools import wraps
import re
import subprocess
import time

from sqlalchemy import desc, distinct
from sqlalchemy.orm import scoped_session

import bottle
from bottle import route, response, request

from sldb.exporting.clone_export import CloneExport
from sldb.exporting.sequence_export import SequenceExport
from sldb.exporting.mutation_export import MutationExporter
import sldb.api.queries as queries
from sldb.common.models import CloneStats, ModificationLog, Sequence
from sldb.exporting.writers import (CLIPWriter, FASTAWriter, FASTQWriter,
                                    CSVWriter)
import sldb.util.lookups as lookups
from sldb.util.nested_writer import NestedCSVWriter


class EnableCors(object):
    """A class to enable Cross-Origin Resource Sharing to facilitate AJAX
    requests.

    """
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


app = bottle.default_app()
app.install(EnableCors())


def with_session(f):
    @wraps(f)
    def _wrapper(*args, **kwargs):
        session = scoped_session(app.config['session_maker'])()
        try:
            return f(session, *args, **kwargs)
        except:
            raise
        finally:
            session.close()

    return _wrapper


def create_response(j=None, code=200, ctype='application/json'):
    if j is None:
        bottle.response.status = code
        return
    bottle.response.content_type = ctype
    bottle.response.status = code
    if ctype == 'application/json':
        return json.dumps(j)
    return j


def get_paging():
    data = bottle.request.json or {}
    return map(int, (
        data.get('page', 1),
        data.get('per_page', 10)
    ))


def decode_run_length(encoding):
    ids = []
    offset = 1
    for match in re.finditer('(T|F)(\d+)', encoding.upper()):
        size = int(match.group(2))
        if match.group(1) == 'T':
            ids.extend(range(offset, offset + size))
        offset += size
    return ids


@app.route('/samples/list', method=['POST', 'OPTIONS'])
@with_session
def samples_list(session):
    return create_response(queries.get_all_studies(session))


@app.route('/sequences/list', method=['POST', 'OPTIONS'])
@with_session
def sequences_list(session):
    fields = bottle.request.json or {}
    return create_response(queries.get_all_sequences(
        session,
        fields.get('filters', {}),
        fields.get('order_field', None),
        fields.get('order_dir', 'asc'),
        get_paging()))


@app.route('/sequence/<sample_id>/<seq_id>', method=['POST', 'OPTIONS'])
@with_session
def sequence(session, sample_id, seq_id):
    return create_response(queries.get_sequence(session, sample_id, seq_id))


@app.route('/clones/list', method=['POST', 'OPTIONS'])
@with_session
def clones_list(session):
    fields = bottle.request.json or {}
    return create_response(queries.get_all_clones(
        session,
        fields.get('filters', {}),
        fields.get('order_field', None),
        fields.get('order_dir', 'asc'),
        get_paging()))


@app.route('/clone/<clone_id>', method=['POST', 'OPTIONS'])
@with_session
def clone(session, clone_id):
    return create_response(queries.get_clone(session, clone_id))


@app.route('/clone/sequences/<clone_id>', method=['POST', 'OPTIONS'])
@with_session
def clone_sequences(session, clone_id):
    fields = bottle.request.json or {}
    return create_response(queries.get_clone_sequences(
        session, clone_id, fields.get('get_collapse', False), get_paging()
    ))


@app.route('/clone/mutations/<clone_id>', method=['POST', 'OPTIONS'])
@with_session
def clone_mutations(session, clone_id):
    fields = bottle.request.json or {}
    return create_response(queries.get_clone_mutations(
        session,
        clone_id,
        fields.get('type', 'percent'),
        int(fields.get('value', 0))
    ))


@app.route('/clone/lineage/<clone_id>', method=['POST', 'OPTIONS'])
@with_session
def clone_lineage(session, clone_id):
    return create_response(queries.get_clone_tree(session, clone_id))


@app.route('/samples/analyze/<sample_encoding>', method=['POST', 'OPTIONS'])
@with_session
def analyze_samples(session, sample_encoding):
    fields = bottle.request.json or {}
    return create_response(
        queries.analyze_samples(
            session,
            decode_run_length(sample_encoding),
            fields.get('filter_type', 'unique_multiple'),
            fields.get('include_outliers', True),
            fields.get('include_partials', True),
            fields.get('percentages', False),
            fields.get('grouping', 'name'),
        )
    )


@app.route('/samples/overlap/<sample_encoding>', method=['POST', 'OPTIONS'])
@with_session
def overlap(session, sample_encoding):
    fields = bottle.request.json or {}
    return create_response(queries.get_clone_overlap(
        session,
        decode_run_length(sample_encoding),
        fields.get('filter_type', 'clones_all'),
        get_paging())
    )


@app.route('/samples/v_usage/<sample_encoding>', method=['POST', 'OPTIONS'])
@with_session
def v_usage(session, sample_encoding):
    fields = bottle.request.json or {}
    data, x_categories, totals = queries.get_v_usage(
        session,
        decode_run_length(sample_encoding),
        fields.get('filter_type', 'unique_multiple'),
        fields.get('include_outliers', True),
        fields.get('include_partials', True),
        fields.get('grouping', 'name'),
        fields.get('by_family', False)
    )

    x_categories.sort()
    y_categories = sorted(data.keys())
    array = []
    min_v = 100
    max_v = 0
    for i, x in enumerate(x_categories):
        for j, y in enumerate(y_categories):
            usage_for_y = data[y]
            if x in usage_for_y:
                min_v = min(min_v, usage_for_y[x])
                max_v = max(max_v, usage_for_y[x])
                array.append([i, j, usage_for_y[x]])
            else:
                array.append([i, j, 0])

    return create_response({
        'x_categories': map(lambda e: 'IGHV{}'.format(e), x_categories),
        'y_categories': y_categories,
        'totals': totals,
        'data': array,
        'min': min_v,
        'max': max_v
    })


def run_rest_service(session_maker, args):
    app.config['session_maker'] = session_maker
    app.run(host='0.0.0.0',
            port=args.port,
            server='gevent',
            debug=args.debug)
