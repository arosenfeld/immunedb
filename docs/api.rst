API --- ``sldb.api``
================================

REST Interface --- ``sldb.api.rest_service``
--------------
The `REST <http://en.wikipedia.org/wiki/Representational_state_transfer>`_
interface is a miniature HTTP server which serves JSON stanzas based on the
requested URL.  For example, requesting ``/api/studies`` will provide a list of
all studies in the database.

The primary use of this is to provide an API for AJAX requests from a
web-interface to the sequence database without the need for writing SQL.

The ``sldb.api.rest_service`` directly handles incoming connections using
`gevent <http://www.gevent.org>`_, parses the URL and then calls the required
functions in :py:class:`sldb.api.queries`.

.. automodule:: sldb.api.rest_service
    :members:
    :exclude-members: relationship

Queries --- ``sldb.api.queries``
-------
.. automodule:: sldb.api.queries
    :members:
    :exclude-members: relationship, desc, distinct

Exporting --- ``sldb.api.export``
---------
The ``sldb.api.export`` module handles exporting of sequences in a variety of
formats, including tab-delimited and FASTA.

This is generally invoked through the REST API, but can also be done manually
with:

.. code-block:: python
    
    from sldb.api.export import SequenceExport

    session = ...
    eformat = 'tab'
    rtype = 'sample'
    rids = [1,2]
    selected_fields = ['seq_id', 'sample_name', 'sequence']
    export = SequenceExport(session, eformat, rtype, rids, selected_fields,
                            duplicates=True, noresults=False)

    for line in export.get_data():
        print line

This will print the sequences from samples 1 and 2 as a tab-delimited file with
fields ``seq_id``, ``sample_name``, and ``sequence`` including duplicate
sequences but excluding those which could not be identified with a V- and
J-gene.

.. automodule:: sldb.api.export
    :members:
