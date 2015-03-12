Overview
============
Terminology
-----------
In this documentation the term **database** will always refer to a mySQL database
and **instance** or **SLDB instance** will always refer to a deployment of SLDB.

This is an important distinction as each SLDB instance utilizes at least two
database: the **master** database and one or more **data** database.

The reason for this is comparability.  Many times it is desirable to analyze the
same data in multiple ways, by either varying the implementation of or parameters
for different steps of the analysis pipeline.  For example, one may want to
identify V and J genes using SLDB's method as well as HighV-Quest's, or assign
clones only to full reads rather than all reads.

In SLDB each of these variations is termed a **version** of the data.  The
master database holds information which does not change across version
[#clone_groups]_s where the data databases contain that which does.
Specifically the master database contains information about:

- **Studies:** Study name, metadata, etc.
- **Samples:** Sample name, date, tissue, etc.
- **Subjects:** Identifier, metadata, etc.

Where the data database contains information about:

- **Sequences:** V, J, CDR3, etc.
- **Clones:** Members, CDR3, etc.
- **Statistics:** Clone frequencies, aggregate sample data, etc.

As more versions of the SLDB instance are added, more data databases will be
created, but the master will remain constant [#clone_groups]_.  This provides
for comparability as IDs and ordering of records in the master remain constant.

It is important to note that new data should be added to all versions of a given
SLDB instance and that entirely different sets of data should be placed in a
separate SLDB instance (with its own master and data databases).


Database Installation
---------------------
SLDB utilizes `mySQL <mysql.com>`_ as its underlying data store.  As it is a
drop-in replacement for mySQL, `MariaDB <mariadb.org>`_ can also be used and is
recommended.  Please see the associated website for installation instructions.

SLDB requires that `TokuDB <tokutek.com/tokudb-for-mysql>`_ be used as the
storage engine for mySQL.  It `has been shown
<http://www.tokutek.com/tokudb-for-mysql/benchmarks-vs-innodb-hdd/>`_ to perform
better on large datasets than the default storage engine InnoDB.  See their
website for installation instructions.

SLDB Installation
-----------------

virtualenv
^^^^^^^^^^

It is recommended that SLDB be installed within a virtual environment, creating
an isolated environment from the rest of the system.

If the `virtualenv package <https://pypi.python.org/pypi/virtualenv>`_ is not
installed, get it with:

.. code-block:: bash

    $ pip install virtualenv

The create and activate a new virtual environment:

.. code-block:: bash

    $ virtualenv sldb-venv
    $ cd sldb-venv
    $ source bin/activate

Finally, get and install SLDB:

.. code-block:: bash

    $ git clone https://github.com/arosenfeld/sldb.git
    $ cd sldb
    $ pip install numpy scipy
    $ python setup.py install

Global
^^^^^^^^^^
.. warning::
    Globally installing SLDB is generally not recommended.  Doing so can result
    in version conflicts and requires root permissions.  Only in specialized
    situations (within a VM) should this be used.

If instead a global install is desired, run:

.. code-block:: bash

    $ git clone https://github.com/arosenfeld/sldb.git
    $ cd sldb
    $ pip install numpy scipy
    $ python setup.py install
