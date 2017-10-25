import json
from functools import wraps
import re
import sys
import time

from sqlalchemy import distinct
from sqlalchemy.orm import scoped_session

import bottle
from bottle import response, request

try:
    from rollbar.contrib.bottle import RollbarBottleReporter
    ROLLBAR_SUPPORT = True
except ImportError:
    ROLLBAR_SUPPORT = False

from immunedb.exporting.clone_export import CloneExport
from immunedb.exporting.sequence_export import SequenceExport
from immunedb.exporting.mutation_export import MutationExporter
import immunedb.api.queries as queries
from immunedb.common.models import CloneStats
from immunedb.exporting.writers import (CLIPWriter, FASTAWriter, FASTQWriter,
                                        CSVWriter)
from immunedb.util.log import logger


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
        except Exception as e:
            raise e
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
    return create_response(queries.get_sequence(session, sample_id, seq_id))


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
    return create_response(queries.get_selection_pressure(session, clone_id))


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


def set_download(exp_type, ext):
    time_str = time.strftime('%Y-%m-%d-%H-%M')
    response.headers['Content-Disposition'] = 'attachment;filename={}'.format(
        '{}_{}.{}'.format(exp_type, time_str, ext))


@app.route('/export/sequences/<from_type>/<encoding>', method=['POST'])
@with_session
def export_sequences(session, from_type, encoding):
    fields = request.forms or {}
    eformat = fields.get('format')

    set_download('sequences', eformat.split('-', 1)[0])

    if eformat == 'csv':
        writer = CSVWriter()
    elif eformat.startswith('fasta'):
        writer = FASTAWriter(replaced_sequences='fill' in eformat)
    elif eformat.startswith('fastq'):
        writer = FASTQWriter(replaced_sequences='fill' in eformat)
    elif eformat.startswith('clip'):
        writer = CLIPWriter(replaced_sequences='fill' in eformat)

    if from_type == 'sample':
        ids = decode_run_length(encoding)
    else:
        ids = [int(encoding)]

    export = SequenceExport(
        session, writer, from_type, ids,
        selected_fields=fields.get('fields').split(','),
        subject_uniques=fields.get('subject_uniques', False),
        only_with_clones=fields.get('only_with_clones', False)
    )
    for line in export.get_data():
        yield line


@app.route('/export/clones/<from_type>/<encoding>', method=['POST'])
@with_session
def export_clones(session, from_type, encoding):
    fields = request.forms or {}

    set_download('clones', 'csv')

    if from_type == 'sample':
        ids = decode_run_length(encoding)
    else:
        ids = [int(encoding)]

    export = CloneExport(
        session, from_type, ids,
        fields.get('fields').split(',')
    )
    for line in export.get_data():
        yield line


@app.route('/export/mutations/<from_type>/<encoding>', method=['POST'])
@with_session
def export_mutations(session, from_type, encoding):
    fields = request.forms or {}

    set_download('mutations', 'csv')

    if from_type == 'sample':
        query = session.query(
            distinct(CloneStats.clone_id).label('clone_id')
        ).filter(
            CloneStats.sample_id.in_(decode_run_length(encoding))
        )
        clone_ids = map(lambda r: r.clone_id, query.all())
    else:
        clone_ids = [int(encoding)]

    export = MutationExporter(
        session, clone_ids,
        decode_run_length(encoding)
        if fields.get('only_sample_rows') else None,
        fields.get('thresh_type', 'sequences'),
        int(fields.get('thresh_value', 0))
    )

    for line in export.get_data():
        yield line


@app.route('/shutdown', method=['POST'])
def shutdown():
    if app.config['allow_shutdown']:
        logger.warning('Shutting down from remote request')
        sys.exit()
    return create_response(code=404)


def run_rest_service(session_maker, args):
    if args.rollbar_token:
        if not ROLLBAR_SUPPORT:
            logger.error('Rollbar is not installed')
            return
        rbr = RollbarBottleReporter(
            access_token=args.rollbar_token,
            environment=args.rollbar_env)
        bottle.install(rbr)

    app.config['session_maker'] = session_maker
    app.config['allow_shutdown'] = args.allow_shutdown
    if args.debug:
        app.catchall = False
    app.run(host='0.0.0.0', port=args.port, server='gevent', debug=args.debug)
