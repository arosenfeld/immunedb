Overview
============
Terminology
-----------
In this documentation the term **database** will always refer to a MySQL database
and **instance** or **SLDB instance** will always refer to a deployment of SLDB.

Dependency Installation
---------------------
MySQL
^^^^^
SLDB utilizes `MySQL <http://mysql.com>`_ as its underlying data store.  As it
is a drop-in replacement for MySQL, `MariaDB <http://mariadb.org>`_ can also be
used and is recommended.  Please see the associated website for installation
instructions.

Baseline & R
^^^^^^^^^^^^
SLDB uses `Baseline <http://selection.med.yale.edu/baseline>`_ to
calculate selection pressure on clones.  This requires `R
<http://www.r-project.org>`_ to be installed along with the `ade4
<http://cran.r-project.org/web/packages/ade4/index.html>`_ package.
Installation is platform dependent.

The newest version of Baseline can be downloaded `here
<http://selection.med.yale.edu/baseline>`_.  The path to the main script will be
needed for clone statistics generation as described in :ref:`stats_generation`.

Clearcut
^^^^^^^^
`Clearcut <http://bioinformatics.hungry.com/clearcut>`_ is used to generate
lineage trees for clones.  After downloading and compiling per the instructions,
note the path to the ``clearcut`` executable which will be required for
generating trees in :ref:`tree_generation`.

Diversity & Haskell
^^^^^^^^^^^^^^^^^^^
SLDB can calculate rarefaction for samples with the `diversity package
<https://hackage.haskell.org/package/diversity>`_ for `Haskell
<https://www.haskell.org>`_.  This is only required if the ``sldb_rest`` command
will be used as detailed in :ref:`supplemental_tools`.


SLDB Installation
-----------------

virtualenv
^^^^^^^^^^

It is recommended that SLDB be installed within a virtual environment, creating
an isolated environment from the rest of the system.

If the `virtualenv package <https://pypi.python.org/pypi/virtualenv>`_ is not
installed, install it with:

.. code-block:: bash

    $ pip install virtualenv

Then create and activate a new virtual environment:

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

.. rubric:: Footnotes

.. [#clone_groups]
    With the exception of the ``clone_groups`` table which will potentially
    change.
