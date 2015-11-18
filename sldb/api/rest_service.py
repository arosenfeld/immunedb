import json
from functools import wraps
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

'''
def check_missing_fields(*fields):
    @wraps(f)
    def _wrapper(*args, **kwargs):
        data = bottle.request.json
        missing = [
            f for f in fields
            if f not in data or (type(data[f]) in (str, unicode) and
                len(data[f].strip()) == 0)
        ]
        if len(missing) > 0:
            return create_response({
                'missing_fields': missing
            }, code=422)
        return f(*args, **kwargs)
    return _wrapper
'''

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
def samples_list(session, sample_id, seq_id):
    return create_response(queries.get_sequence(session, sample_id, seq_id))


def run_rest_service(session_maker, args):
    app.config['session_maker'] = session_maker
    app.run(host='0.0.0.0',
            port=args.port,
            server='gevent',
            debug=args.debug)
