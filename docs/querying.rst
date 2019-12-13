.. _querying:

Querying with SQL
=================
ImmuneDB is backed by a MySQL database that can be queried directly to gather
information, bypassing the Python API.

Accessing the Database
----------------------
There are many ways to access the database directly.  The two introduced here
are directly through MySQL or using ``immunedb_sql`` which simply wraps a call to
MySQL.

With the ``immunedb_sql`` wrapper (recommended)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    $ immunedb_sql PATH_TO_CONFIG

This is entirely equivalent to using ``mysql`` and will drop to the MySQL
interpreter.  You can also pass a query directly from the command line.  For
example:

.. code-block:: bash

    $ immunedb_sql PATH_TO_CONFIG --query 'select * from samples'


Directly with MySQL
^^^^^^^^^^^^^^^^^^^

From the command line, you may access an ImmuneDB database ``DATABASE`` from user
``USERNAME`` with:

.. code-block:: bash

    $ mysql -u USERNAME -p DATABASE

This will prompt for a password and then to the database.  This method of access
is useful for quickly querying the database.  To save results of a query
``QUERY`` run the command:

.. code-block:: bash

    $ mysql -u USERNAME -p DATABASE -e "QUERY" > output
