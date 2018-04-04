Directly Querying the Database
==============================
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

Since ImmuneDB stores usernames and passwords in config files ``immunedb_sql`` is provided
as a small wrapper around the ``mysql`` command.  It can be invoked with:


Querying
--------

The data is split into different referential tables based on the data models
defined in :ref:`Models <models>`.

Most tables contain data that will generally be aggregated by a query directly.
For example, to determine how many sequences are in sample 10, one could simply
count the rows:

.. code-block:: sql

    > SELECT COUNT(*) FROM sequences WHERE sample_id=10;

This query will work, but may be slow with large datasets.  Two tables can
assist in some computationally-expensive queries by providing pre-aggregated
information.  The ``CloneStats`` and ``SampleStats`` tables are pre-populated
with the ``immunedb_clone_stats`` and ``immunedb_sample_stats`` :doc:`pipeline
</pipeline>` commands.

This same result could be achieved more quickly with:

.. code-block:: sql

    > SELECT sequence_cnt FROM sample_stats WHERE filter_type='all' AND
        outlier=1 AND full_reads=0 AND sample_id=10;

This is more verbose but is quicker and also makes other more complex tasks
simpler.  For example, let's find how many unique sequences are in sample 10
which are full reads and not outliers:

.. code-block:: sql

    > SELECT sequence_cnt FROM sample_stats WHERE filter_type='unique' AND
        outlier=0 AND full_reads=1 AND sample_id=10;

The ``CloneStats`` table also has useful information that has been
pre-aggregated.  For example, how many unique sequences in sample 10 are in
clone 5:

.. code-block:: sql

    > SELECT unique_cnt FROM clone_stats WHERE sample_id=10 AND clone_id=5;

Or how many total unique sequences are in clone 5 (``NULL`` is a placeholder in
the ``sample_id`` column meaning "All Samples"):

.. code-block:: sql

    > SELECT unique_cnt FROM clone_stats WHERE sample_id=NULL AND clone_id=5;


Other Example Queries
---------------

**How many clones have a CDR3 starting with `CARD`?**

.. code-block:: sql

    > SELECT COUNT(*) FROM clones WHERE clones.cdr3_aa LIKE 'CARD%';

**Get a list of non-identifiable sequences in FASTA format.**

.. code-block:: sql

    > SELECT CONCAT('>', seq_id, '\n', sequence) FROM noresults;

**How many indels and total sequences do I have?"**

.. code-block:: sql

    > SELECT SUM(IF(probable_indel_or_misalign=1, 1, 0)) AS indels, COUNT(*) AS
    total FROM sequences;
