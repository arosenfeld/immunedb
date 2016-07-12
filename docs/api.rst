API --- ``immunedb.api``
====================

.. note::

    This documentation is still in-the-works and is incomplete.

REST Interface --- ``immunedb.api.rest_service``
--------------------------------------------
The `REST <http://en.wikipedia.org/wiki/Representational_state_transfer>`_
interface is a miniature HTTP server which serves JSON stanzas based on the
requested URL.  For example, requesting ``/api/studies`` will provide a list of
all studies in the database.

The primary use of this is to provide an API for AJAX requests from a
web-interface to the sequence database without the need for writing SQL.

The ``immunedb.api.rest_service`` directly handles incoming connections using
`gevent <http://www.gevent.org>`_, parses the URL and then calls the required
functions in :py:class:`immunedb.api.queries`.

.. automodule:: immunedb.api.rest_service
    :members:
    :exclude-members: relationship

Queries --- ``immunedb.api.queries``
--------------------------------
.. automodule:: immunedb.api.queries
    :members:
    :exclude-members: relationship, desc, distinct
