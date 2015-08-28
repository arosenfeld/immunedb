Directly Querying the Database
==============================
SLDB is backed by a MySQL database that can be queried directly to gather
information, bypassing the Python API.

Accessing the Database
----------------------
There are many ways to access the database directly for example through a DBAPI
available for all commonly used languages.  Alternatively, one can access MySQL
directly from the command line with:

.. code-block:: bash

    $ myqsql -u USERNAME -p DATABASE

This will connect to the database ``DATABASE`` with the username ``USERNAME``.

This method of access is useful for quickly querying the database.  To save
results of a query ``QUERY`` run the command:

.. code-block:: bash

    $ myqsql -u USERNAME -p DATABASE -e "QUERY" > output


Querying
--------

The data is split into different referential tables based on the data models
defined in :doc:`Models </models>`.

Most tables contain data that will generally be aggregated by a query directly.
For example, to determine how many sequences are in sample 10, one could simply
count the rows:

.. code-block:: sql

    > SELECT COUNT(*) FROM sequences where sample_id=10;

This query will work, but may be slow with large datasets.  Two tables can
assist in some computationally-expensive queries by providing pre-aggregated
information.  The ``CloneStats`` and ``SampleStats`` tables are pre-populated
with the ``sldb_clone_stats`` and ``sldb_sample_stats`` :doc:`pipeline
</pipeline>` commands.

This same result could be achieved more quickly with:

.. code-block:: sql

    > SELECT sequence_cnt FROM sample_stats WHERE filter_type='all' AND
        outlier=1 AND full_reads=0 and sample_id=10;

This is more verbose but is quicker and also makes other more complex tasks
simpler.  For example, let's find how many unique sequences are in sample 10
which are full reads and not outliers:

.. code-block:: sql

    > SELECT sequence_cnt FROM sample_stats WHERE filter_type='unique' AND
        outlier=0 AND full_reads=1 and sample_id=10;

The ``CloneStats`` table also has useful information that has been
pre-aggregated.  For example, how many unique sequences in sample 10 are in
clone 5:

.. code-block:: sql

    > SELECT unique_cnt FROM clone_stats where sample_id=10 and clone_id=5;

Or how many total unique sequences are in clone 5 (``NULL`` is a placeholder in
the ``sample_id`` column meaning "All Samples"):

.. code-block:: sql

    > SELECT unique_cnt FROM clone_stats where sample_id=NULL and clone_id=5;


Other Example Queries
---------------

**How many clones have a CDR3 starting with `CARD`?**

.. code-block:: sql

    > SELECT COUNT(*) FROM clones WHERE clones.cdr3_aa like 'CARD%';

**Get a list of non-identifiable sequences in FASTA format.**

.. code-block:: sql

    > SELECT CONCAT('>', seq_id, '\n', sequence) from noresults;

.. note::
    To output to a file try ``mysql -r -u USER -p DATABASE -e "`SELECT
    CONCAT('>', seq_id, '\n', sequence) from noresults`" > output.fasta``

**How many indels and total sequences do I have?"**

.. code-block:: sql

    > SELECT SUM(IF(probable_indel_or_misalign=1, 1, 0)) AS indels, count(*) AS
    total FROM sequences;
