import argparse
import json
import math
import time

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, scoped_session

import bottle
from bottle import route, response, request, install, run

from sldb.api.export import CloneExport, SequenceExport
import sldb.api.queries as queries
from sldb.common.models import *


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
    """Gets a list of all sequences.

    :returns: A list of all sequences
    :rtype: str

    """
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
    """Gets the sequence identified by ``seq_id`` in sample with id
    ``sample_id``.

    :param int sample_id: The sample ID of the sequence
    :param str seq_id: The sequence ID of the sequence

    :returns: The requested sequence if it exists
    :rtype: str

    """
    session = scoped_session(session_factory)()
    seq = queries.get_sequence(session, int(sample_id), seq_id)
    session.close()
    return json.dumps({'sequence': seq})


@route('/api/studies')
def studies():
    """Gets a list of all studies and their associated samples.

    :returns: A list of all studies and their associated samples
    :rtype: str

    """
    session = scoped_session(session_factory)()
    studies = queries.get_all_studies(session)
    session.close()
    return json.dumps({'studies': studies})


@route('/api/subjects')
def subjects():
    """Gets a list of all subjects.

    :returns: A list of all subjects
    :rtype: str

    """
    session = scoped_session(session_factory)()
    subjects = queries.get_all_subjects(session, _get_paging())
    session.close()
    return json.dumps({'subjects': subjects})


@route('/api/subject/<sid>')
def subject(sid):
    """Gets the subject with id ``sid``.

    :param int sid: The subject ID to query

    :returns: The requested subject if it exists
    :rtype: str

    """
    session = scoped_session(session_factory)()
    subject = queries.get_subject(session, int(sid))
    session.close()
    return json.dumps({'subject': subject})


@route('/api/clones/')
def clones():
    """Gets a list of all clones.

    :returns: A list of all clones
    :rtype: str

    """
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
    """Compares clones, outputting their mutations and other pertinent
    information.

    :param str uids: A list of clones to compare either as just their IDs or \
    ``ID_SAMPLE`` to limit that clone to a given sample

    :returns: A clone comparison between ``uids``
    :rtype: str

    """
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


@route('/api/clone_tree/<cid>')
def clone_tree(cid):
    """ Gets the lineage tree represented by JSON for a clone.

    :param int cid: The clone ID of which to get the lineage tree

    :returns: The lineage tree as JSON
    :rtype: str

    """
    session = scoped_session(session_factory)()
    tree = queries.get_clone_tree(session, cid)
    session.close()
    return tree


@route('/api/clone_overlap/<filter_type>/<samples>')
@route('/api/subject_clones/<filter_type>/<subject>')
def clone_overlap(filter_type, samples=None, subject=None):
    """Gets the clones that overlap between a set of samples. If ``samples`` is
    supplied, the overlap of clones between those is returned.  If ``samples``
    is not supplied, the clonal overlap for all samples from ``subject`` is
    returned.

    :param str filter_type: The filter for the clones to apply.  This can \
    currently be ``clones_all`` for all clones, ``clones_functional`` for \
    only    functional clones, or ``clones_nonfunctional`` for only \
    non-functional clones.
    :param str samples: A comma-separated sequences of sample IDs if \
    comparing specific samples.  Otherwise, ``None``.
    :param int subject: The subject ID to use for clonal overlap.  If \
    specified, all samples from the subject are compared.

    :returns: The overlap of clones in the specified sample or subjects
    :rtype: str

    """
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
    """Gets the statistics for a given set of samples both including and
    excluding outliers.

    :param str samples: A comma-separated list of sample IDs for which to \
    gather statistics

    :returns: Statistics for all samples in ``samples``
    :rtype: str

    """
    session = scoped_session(session_factory)()
    samples = _split(samples)
    ret = json.dumps({
        'outliers': {
            'all_reads': queries.get_stats(session, samples,
                                           include_outliers=True,
                                           full_reads=False),
            'full_reads': queries.get_stats(session, samples,
                                            include_outliers=True,
                                            full_reads=True),
        },
        'no_outliers': {
            'all_reads': queries.get_stats(session, samples,
                                           include_outliers=False,
                                           full_reads=False),
            'full_reads': queries.get_stats(session, samples,
                                            include_outliers=False,
                                            full_reads=True),
        }
    })
    session.close()
    return ret


@route('/api/modification_log')
def modification_log():
    session = scoped_session(session_factory)()
    logs = []
    for log in session.query(ModificationLog).order_by(
            ModificationLog.datetime):
        logs.append({
            'datetime': log.datetime.strftime('%Y-%m-%d %H:%M:%S'),
            'action_type': log.action_type,
            'info': json.loads(log.info),
        })
    session.close()
    return json.dumps({'logs': logs})


@route('/api/v_usage/<filter_type>/<outliers>/<full_reads>/<samples>')
def v_usage(filter_type, outliers, full_reads, samples):
    """Gets the V usage for samples in a heatmap-formatted array.

    :param str filter_type: The filter type of sequences for the v_usage
    :param str samples: A comma-separated string of sample IDs for v_usage

    :returns: The V usage as an [(x,y)...] JSON array along with x and y
    categories
    :rtype: str

    """
    session = scoped_session(session_factory)()
    data, headers = queries.get_v_usage(
        session, _split(samples), filter_type,
        outliers == 'true', full_reads == 'true')
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


@route('/api/data/export_clones/<rtype>/<rids>', methods=['GET'])
@route('/api/data/export_clones/<rtype>/<rids>/', methods=['GET'])
def export_clones(rtype, rids):
    """Downloads a tab-delimited file of clones.

    :param str rtype: The type of record to filter the query on.  Currently
        either "sample" or "clone"
    :param str rids: A comma-separated list of IDs of ``rtype`` to export

    :returns: The properly formatted export data
    :rtype: str

    """
    assert rtype in ('sample', 'clone')

    session = scoped_session(session_factory)()
    fields = _get_arg('fields', False).split(',')
    include_total_row = _get_arg('include_total_row', False) == 'true' or False

    name = '{}_{}.tab'.format(
        rtype,
        time.strftime('%Y-%m-%d-%H-%M'))

    response.headers['Content-Disposition'] = 'attachment;filename={}'.format(
        name)

    export = CloneExport(session, rtype, _split(rids), fields,
                         include_total_row)
    for line in export.get_data():
        yield line

    session.close()


@route('/api/data/export_sequences/<eformat>/<rtype>/<rids>', methods=['GET'])
@route('/api/data/export_sequences/<eformat>/<rtype>/<rids>/', methods=['GET'])
def export_sequences(eformat, rtype, rids):
    """Downloads an exported format of specified sequences.

    :param str eformat: The export format to use.  Currently "tab", "orig", and
        "clip" for tab-delimited, FASTA, FASTA with filled in germlines, and
        FASTA in CLIP format respectively
    :param str rtype: The type of record to filter the query on.  Currently
        either "sample" or "clone"
    :param str rids: A comma-separated list of IDs of ``rtype`` to export

    :returns: The properly formatted export data
    :rtype: str

    """

    assert eformat in ('tab', 'fill', 'orig', 'clip')
    assert rtype in ('sample', 'clone')

    session = scoped_session(session_factory)()

    fields = _get_arg('fields', False).split(',')
    min_copy = _get_arg('min_copy_number', False)
    min_copy = int(min_copy) if min_copy is not None else 1

    if eformat == 'tab':
        name = '{}_{}.tab'.format(
            rtype,
            time.strftime('%Y-%m-%d-%H-%M'))
    else:
        if 'seq_id' in fields:
            fields.remove('seq_id')
        name = '{}_{}_{}.fasta'.format(
            rtype,
            eformat,
            time.strftime('%Y-%m-%d-%H-%M'))

    response.headers['Content-Disposition'] = 'attachment;filename={}'.format(
        name)

    export = SequenceExport(
        session, eformat, rtype, _split(rids), fields,
        min_copy=min_copy,
        duplicates=_get_arg('duplicates', False) == 'true',
        noresults=_get_arg('noresults', False) == 'true')
    for line in export.get_data():
        yield line

    session.close()


def run_rest_service(session_maker, args):
    """Runs the rest service based on command line arguments"""
    global session_factory
    session_factory = session_maker
    bottle.install(EnableCors())
    if args.debug:
        bottle.run(host='0.0.0.0', port=args.port, server='gevent',
                   debug=True)
    else:
        bottle.run(host='0.0.0.0', port=args.port, server='gevent')
