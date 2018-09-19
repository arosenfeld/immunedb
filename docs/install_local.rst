Installing Locally (advanced)
*****************************
This section details how to set ImmuneDB up locally on a machine.  This is a
more complicated process than using the :doc:`Docker method <install_docker>`
but may be useful if you plan on running ImmuneDB remotely on a server rather
than locally.

Dependency Installation
=======================

MySQL
-----
ImmuneDB utilizes `MySQL <http://mysql.com>`_ as its underlying data store.  We
recommend using its drop-in replacement, `MariaDB <http://mariadb.org>`_.
Please consult their website and your operating systems package manager for
installation instructions.

R (optional)
------------
`Baseline <http://selection.med.yale.edu/baseline>`_ can optionally be used to
calculate selection pressure on clones.  This requires `R
<http://www.r-project.org>`_ to be installed along with the `ade4
<http://cran.r-project.org/web/pack:ges/ade4/index.html>`_ package.
Installation is platform dependent.

The newest version of Baseline can be downloaded `here
<http://selection.med.yale.edu/baseline>`_.  The path to the main script will
be needed for clone statistics generation as described in
:ref:`stats_generation`.

For genotyping, `TIgGER <http://tigger.readthedocs.io>`_ must also be
installed.

Bowtie2 (optional)
------------------
`Bowtie2 <bowtie-bio.sourceforge.net>`_ can be used to locally align sequences
which cannot be aligned using the built-in anchor method.

Clearcut (optional)
-------------------
`Clearcut <http://bioinformatics.hungry.com/clearcut>`_ can be used to generate
lineage trees for clones.  After downloading and compiling per the instructions,
note the path to the ``clearcut`` executable which will be required for
generating trees in :ref:`tree_generation`.

ImmuneDB Installation
=====================

It is recommended that ImmuneDB be installed within a `venv`, creating
an isolated environment from the rest of the system.

To create a virtual environment and activate it run:

.. code-block:: bash

    $ python3 -m venv immunedb
    $ source immunedb/bin/activate

Then install ImmuneDB:

.. code-block:: bash

    $ pip install immunedb

Web Interface Installation
==========================
Please refer to the `ImmuneDB Frontend installation instructions
<https://github.com/arosenfeld/immunedb-frontend#immunedb-frontend>`_.
