.. _api:

Python API
==========

.. note::

    This section is currently incomplete.  We're working to fill out the
    details of the Python API as soon as possible.

Configuration
-------------
The ``immunedb.common.config`` module provides methods to initialize a
connection to a new or existing database.

Most programs using ImmuneDB will start with code similar to:

.. code-block:: python


    import immunedb.common.config as config


    parser = config.get_base_arg_parser('Some description of the program')
    # ... add any additional arguments to the parser ...
    args = parser.parse_args()

    session = config.init_db(args.db_config)

When this script is run, it will require at least one argument which is the
path to a database configuration (as generated with ``immunedb_admin``).  Using
that, a ``Session`` object will be made, connected to the associated database.

One can also directly specify the path to a configuration directly.

.. code-block:: python


    import immunedb.common.config as config


    session = config.init_db('path/to/config')

Alternatively a dictionary with the same information can be passed:

.. code-block:: python


    import immunedb.common.config as config


    session = config.init_db({
        'host': '...',
        'database': '...',
        'username': '...',
        'password': '...',
    })

Returned will be a ``Session`` object which can be used to interact with the
database.

Using the Session
-----------------
ImmuneDB is built using `SQLAlchemy <http://sqlalchemy.org>`_ as a MySQL
abstraction layer.  Simply put, instead of writing SQL, the database is queried
using Python constructs.  Full documentation on using the session can be found
in `SQLAlchemy's documentation
<http://docs.sqlalchemy.org/en/latest/orm/session.html>`_.

Once a session is created, the models listed below can be queried.

Example Queries
---------------
Below are some example queries that demonstrate how to use the ImmuneDB API.

Clone CDR3s
^^^^^^^^^^^
Get all clones with a given V-gene and print their CDR3 AA
sequences.

**Input**

.. code-block:: python

    import immunedb.common.config as config
    from immunedb.common.models import Clone

    session = config.init_db(...)

    for clone in session.query(Clone).filter(Clone.v_gene == 'IGHV3-30'):
        print('clone {} has AAs {}'.format(clone.id, clone.cdr3_aa))

**Output**

.. code-block:: bash

    clone 37884 has AAs CARGYSSSYFDYW
    clone 37886 has AAs CARSRTSLSIYGVVPTGDFDSW
    clone 37885 has AAs CARNGLNTVSGVVISPKYWLDPW
    clone 37887 has AAs CARDLFRGVDFYYYGMDVW

Clone Frequency
^^^^^^^^^^^^^^^
Determine how many sequences appear in each sample belonging to clone 1234.

Note the ``CloneStats`` model has one entry for each clone/sample combination
plus one where the ``sample_id`` field is ``null`` which represents the overall
clone.

**Input**

.. code-block:: python

    import immunedb.common.config as config
    from immunedb.common.models import CloneStats

    session = config.init_db(...)
    for stat in session.query(CloneStats).filter(
            CloneStats.clone_id == 1234).order_by(CloneStats.sample_id):
        print('clone {} has {} unique sequences and {} copies {}'.format(
            stat.clone_id,
            stat.unique_cnt,
            stat.total_cnt,
            ('in sample ' + stat.sample.name) if stat.sample else 'overall'))

**Output**

.. code-block:: bash

    clone 1234 has 53 unique sequences and 1331 copies overall
    clone 1234 has 27 unique sequences and 379 copies in sample sample1
    clone 1234 has 27 unique sequences and 339 copies in sample sample3
    clone 1234 has 24 unique sequences and 311 copies in sample sample4
    clone 1234 has 28 unique sequences and 302 copies in sample sample10

V-gene Usage
^^^^^^^^^^^^
This is a more complex query which gathers the V-gene usage of all sequences
which are (a) in subject with ID 1, (b) associated with a clone, and (c) are
unique to the subject, printing them from least to most frequent.

**Input**

.. code-block:: python

    import immunedb.common.config as config
    from immunedb.common.models import Sequence, SequenceCollapse

    session = config.init_db(...)

    subject_unique_seqs = session.query(
        func.count(Sequence.seq_id).label('count'),
        Sequence.v_gene
    ).join(
        SequenceCollapse
    ).filter(
        Sequence.subject_id == 1,
        ~Sequence.clone_id.is_(None),
        SequenceCollapse.copy_number_in_subject > 0
    ).group_by(
        Sequence.v_gene
    ).order_by(
        'count'
    )

    for seq in subject_unique_seqs:
        print(seq.v_gene, seq.count)

**Output**

.. code-block:: bash

    # ... output truncated ...
    IGHV4-34 1128
    IGHV1-2 1160
    IGHV3-48 1169
    IGHV4-39 1310
    IGHV3-7 1345
    IGHV3-30|3-30-5|3-33 1607
    IGHV3-23|3-23D 1626
    IGHV3-21 1878

.. _models:

Data Models
-----------

.. automodule:: immunedb.common.models
    :members:
    :exclude-members: relationship
