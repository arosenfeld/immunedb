Getting Started (Bare Metal)
====================================
This section details how to set SLDB up locally on a machine.  Alternatively,
you can use :ref:`Docker Compose <docker_install>` to get up and running more
quickly.

Dependency Installation
---------------------
MySQL
^^^^^
SLDB utilizes `MySQL <http://mysql.com>`_ as its underlying data store.  As it
is a drop-in replacement for MySQL, `MariaDB <http://mariadb.org>`_ can also be
used and is recommended.  Please see the associated website for installation
instructions.

Baseline & R (optional)
^^^^^^^^^^^^^^^^^^^^^^^
SLDB can use `Baseline <http://selection.med.yale.edu/baseline>`_ to calculate
selection pressure on clones.  This requires `R <http://www.r-project.org>`_ to
be installed along with the `ade4
<http://cran.r-project.org/web/pack:ges/ade4/index.html>`_ package.
Installation is platform dependent.

The newest version of Baseline can be downloaded `here
<http://selection.med.yale.edu/baseline>`_.  The path to the main script will be
needed for clone statistics generation as described in :ref:`stats_generation`.

Clearcut (optional)
^^^^^^^^^^^^^^^^^^^
`Clearcut <http://bioinformatics.hungry.com/clearcut>`_ can be used to generate
lineage trees for clones.  After downloading and compiling per the instructions,
note the path to the ``clearcut`` executable which will be required for
generating trees in :ref:`tree_generation`.


SLDB Installation
-----------------

virtualenv
^^^^^^^^^^

It is recommended that SLDB be installed within a virtual environment, creating
an isolated environment from the rest of the system.

If `virtualenv package <https://pypi.python.org/pypi/virtualenv>`_ is not
installed, install it globally with:

.. code-block:: bash

    $ sudo pip install virtualenv

Then create and activate a new virtual environment where `NAME` should be
replaced with an appropriate name:

.. code-block:: bash

    $ virtualenv NAME
    $ cd NAME
    $ source bin/activate

Finally, get and install SLDB:

.. code-block:: bash

    $ git clone https://github.com/arosenfeld/sldb.git
    $ cd sldb
    $ pip install numpy
    $ pip install scipy
    $ python setup.py install

Global
^^^^^^
.. warning::
    Installing many packages globally is not recommended.  Using virtual
    environments keeps dependencies separated from the root filesystem.  Only in
    specialized situations (e.g. within a VM) should SLDB installed globally.

If instead a global install is desired, run:

.. code-block:: bash

    $ git clone https://github.com/arosenfeld/sldb.git
    $ cd sldb
    $ pip install numpy scipy
    $ python setup.py install
