Common Classes & Functions --- ``immunedb.common``
==============================================

Configuration --- ``immunedb.common.config``
----------------------------------------
The ``immunedb.common.config`` module provides methods to initialize a
connection to a new or existing database.

Most programs using ImmuneDB will start with code similar to:

.. code-block:: python


    import immunedb.common.config as config
    parser = config.get_base_arg_parser('Some description of the program')
    # ... add any additional arguments to the parser ...
    args = parser.parse_args()

    session = config.init_db(args.db_config)

This takes the first two command line arguments as the master and data
configuration file paths and initializes the database.  Of course one can also
manually specify the paths with simply:

.. code-block:: python


    import immunedb.common.config as config
    session = config.init_db('path/to/config')

To avoid using a configuration file, the ``from_dict`` argument can be set to
``True``, and the configurations will be read from dictionaries:

.. code-block:: python


    import immunedb.common.config as config
    session = config.init_db({
        'host': '...',
        'database': '...',
        'username': '...',
        'password': '...',
    })

All will return a ``Session`` object which can be used to interact with the
database.

.. automodule:: immunedb.common.config
    :members:

Mutations --- ``immunedb.common.mutations``
---------------------------------------
Mutation statistics are calculated with the :py:class:`Mutations
<immunedb.common.mutations.ContextualMutations>` class.  This is generally done for
clones only but can be done with any set of equal-length sequences for which a
base (germline) sequence can be specified.

.. automodule:: immunedb.common.mutations
    :members:
