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

This takes the first command line arguments as database configuration file path
and initializes the database.  Of course one can also manually specify the
paths with simply:

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

Data Models
-----------

.. automodule:: immunedb.common.models
    :members:
    :exclude-members: relationship
