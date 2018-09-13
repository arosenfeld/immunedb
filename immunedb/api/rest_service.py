import json
from functools import partial, wraps
import os
import re
import signal
import time

from sqlalchemy.orm import scoped_session

import bottle
from bottle import response, request

try:
    from rollbar.contrib.bottle import RollbarBottleReporter
    ROLLBAR_SUPPORT = True
except ImportError:
    ROLLBAR_SUPPORT = False

import immunedb.api.queries as queries
from immunedb.common.models import Sequence, SequenceCollapse
from immunedb.util.log import logger
import immunedb.exporting as exporting


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
    if data.get('nopage', False):
        return None
    return [
        int(data.get('page', 1)),
        int(data.get('per_page', 10))
    ]


def decode_run_length(encoding):
    ids = []
    offset = 1
    for match in re.finditer('(T|F)(\d+)', encoding.upper()):
        size = int(match.group(2))
        if match.group(1) == 'T':
            ids.extend(range(offset, offset + size))
        offset += size
    return ids


def create_app(session_maker, allow_shutdown=False):
    app = bottle.Bottle()
    app.install(EnableCors())

    def with_session(f):
        @wraps(f)
        def _wrapper(*args, **kwargs):
            session = scoped_session(session_maker)()
            try:
                return f(session, *args, **kwargs)
            except Exception:
                raise
            finally:
                session.close()

        return _wrapper

    @app.route('/samples/list', method=['POST', 'OPTIONS'])
    @with_session
    def samples_list(session):
        return create_response(queries.get_samples(session))

    @app.route('/sequences/list', method=['POST', 'OPTIONS'])
    @with_session
    def sequences_list(session):
        fields = bottle.request.json or {}
        return create_response(queries.get_sequences(
            session,
            fields.get('filters', {}),
            fields.get('order_field', None),
            fields.get('order_dir', 'asc'),
            fields.get('subject_id', None),
            get_paging()))

    @app.route('/sequence/<sample_id>/<seq_id>', method=['POST', 'OPTIONS'])
    @with_session
    def sequence(session, sample_id, seq_id):
        return create_response(queries.get_sequence(session, sample_id,
                                                    seq_id))

    @app.route('/clones/list', method=['POST', 'OPTIONS'])
    @app.route('/clones/list/<subject_id>', method=['POST', 'OPTIONS'])
    @with_session
    def clones_list(session, subject_id=None):
        fields = bottle.request.json or {}
        return create_response(queries.get_clones(
            session,
            fields.get('filters', {}),
            fields.get('order_field', None),
            fields.get('order_dir', 'asc'),
            subject_id,
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

    @app.route('/clone/pressure/<clone_id>', method=['POST', 'OPTIONS'])
    @with_session
    def clone_pressure(session, clone_id):
        return create_response(list(
            queries.get_selection_pressure(session, clone_id)))

    @app.route('/samples/analyze/<sample_encoding>', method=[
        'POST', 'OPTIONS'])
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

    @app.route('/samples/overlap/<sample_encoding>', method=[
        'POST', 'OPTIONS'])
    @with_session
    def overlap(session, sample_encoding):
        fields = bottle.request.json or {}
        return create_response(queries.get_clone_overlap(
            session,
            decode_run_length(sample_encoding),
            fields.get('filter_type', 'clones_all'),
            get_paging())
        )

    @app.route('/samples/v_usage/<sample_encoding>', method=[
        'POST', 'OPTIONS'])
    @with_session
    def v_usage(session, sample_encoding):
        fields = bottle.request.json or {}
        data, x_categories, totals, prefix = queries.get_v_usage(
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
            'x_categories': ['{}{}'.format(prefix, e) for e in x_categories],
            'y_categories': y_categories,
            'totals': totals,
            'data': array,
            'min': min_v,
            'max': max_v
        })

    @app.route('/subjects/list', method=['POST', 'OPTIONS'])
    @with_session
    def subjects_list(session):
        return create_response(queries.get_all_subjects(session, get_paging()))

    @app.route('/subject/<subject_id>', method=['POST', 'OPTIONS'])
    @with_session
    def subject(session, subject_id):
        return create_response(queries.get_subject(session, subject_id))

    def set_file(suffix, ext='tsv'):
        time_str = time.strftime('%Y-%m-%d-%H-%M')
        fn = '{}_{}.{}'.format(time_str, suffix, ext)
        response.headers['Content-Disposition'] = 'attachment;filename=' + fn

    @app.route('/export/sequences/<schema>', method=['GET', 'OPTIONS'])
    @with_session
    def export_sequences(session, schema):
        if schema not in exporting.mappings:
            return create_response(code=400)

        seqs = session.query(Sequence)

        min_copies = request.query.get('min_subject_copies', 0)
        if min_copies:
            seqs = seqs.filter(
                SequenceCollapse.copy_number_in_subject >= min_copies
            )

        if request.query.get('clones_only', False):
            seqs = seqs.filter(~Sequence.clone_id.is_(None))

        set_file(schema)

        yield from exporting.get_sequences(session, seqs, schema)

    @app.route('/export/clones/<schema>', method=['GET', 'OPTIONS'])
    @with_session
    def export_clones(session, schema):
        funcs = {
            'summary': partial(
                exporting.get_clone_summary,
                include_lineages=request.query.get('lineages', False)
            ),
            'overlap': exporting.get_clone_overlap,
            'selection': exporting.get_selection,
            'vdjtools': exporting.get_vjdtools_as_zip,
        }
        if schema not in funcs:
            return create_response(code=400)

        set_file('clone_' + schema,
                 'zip' if schema == 'vdjtools' else 'tsv')

        yield from funcs[schema](session)

    @app.route('/export/samples', method=['GET', 'OPTIONS'])
    @with_session
    def export_samples(session):
        set_file('samples')
        yield from exporting.get_samples(session)

    @app.route('/shutdown', method=['GET'])
    def shutdown():
        if allow_shutdown:
            logger.warning('Shutting down from remote request')
            os.kill(os.getppid(), signal.SIGINT)
        return create_response(code=404)

    return app


def run_rest_service(session_maker, args):
    if args.rollbar_token:
        if not ROLLBAR_SUPPORT:
            logger.error('Rollbar is not installed')
            return
        rbr = RollbarBottleReporter(
            access_token=args.rollbar_token,
            environment=args.rollbar_env)
        bottle.install(rbr)

    app = create_app(session_maker, args.allow_shutdown)
    if args.debug:
        app.catchall = False
    app.run(host='0.0.0.0', port=args.port, server=args.server,
            debug=args.debug)
