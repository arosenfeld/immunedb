Common Classes & Functions --- ``sldb.common``
================================

Configuration --- ``sldb.common.config``
-------------
The ``sldb.common.config`` module provides methods to initialize a
connection to a new or existing database.

Most programs using SLDB will start with code similar to:

.. code-block:: python


    import sldb.common.config as config
    parser = config.get_base_arg_parser('Some description of the program')
    # ... add any additional arguments to the parser ...
    args = parser.parse_args()

    session = config.init_db(args.master_db_config, args.data_db_config)

This takes the first two command line arguments as the master and data
configuration file paths and initializes the database.  Of course one can also
manually specify the paths with simply:

.. code-block:: python


    import sldb.common.config as config
    session = config.init_db('path/to/master/config', 'path/to/data/config')

To avoid using a configuration file, the ``from_dict`` argument can be set to
``True``, and the configurations will be read from dictionaries:

.. code-block:: python


    import sldb.common.config as config
    session = config.init_db({
        'host': '...',
        'database': '...',
        'usertname': '...',
        'password': '...',
    }, {
        'host': '...',
        'database': '...',
        'usertname': '...',
        'password': '...',
    }, from_dict=True)

All will return a ``Session`` object which can be used to interact with the
database.

.. WARNING::
    ``config.init_db`` must be called prior to importing any SLDB models from
    the ``sldb.common.models`` package.

.. automodule:: sldb.common.config
    :members:

Mutations --- ``sldb.common.mutations``
---------
Mutation statistics are calculated with the :py:class:`Mutations
<sldb.common.mutations.Mutations>` class.  This is generally done for clones
only but can be done with any set of equal-length sequences for which a base
(germline) sequence can be specified.

.. autoclass:: sldb.common.mutations.Mutations
    :members:

.. autoclass:: sldb.common.mutations.MutationType
    :members:
