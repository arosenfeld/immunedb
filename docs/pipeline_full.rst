Running the Pipeline on Your Data
*********************************

This page describes how to run the ImmuneDB pipeline on raw FASTA/FASTQ data.
It is assumed that you've previously tried the :ref:`example pipeline
<pipeline_example>` and understand the basics of running commands in the Docker
container.

Like in the example, each code block has a header saying if the command should
be run on the host or in the Docker container.

Copying Your Sequence Data Into Docker
======================================
Unlike in the :ref:`example pipeline <pipeline_example>` where sequencing data
was provided, you'll need to copy your own FASTA/FASTQ sequencing data into the
Docker container.

To do so, on the host, we create a new directory in the shared directory
into which we'll copy your sequencing data.  Here we're calling it
``sequences`` but you'll probably want to choose a more descriptive name:

.. code-block:: bash
   :caption: Run on Host

   $ mkdir $HOME/immunedb_share/sequences
   $ cp PATH_TO_SEQUENCES $HOME/immunedb_share/sequences


Creating a Metadata Sheet
==================================
Next, we'll use the ``immunedb_metadata`` command to create a template metadata
file for your sequencing data.  In the Docker container run:

.. code-block:: bash
   :caption: Run in Docker

   $ cd /share/sequences
   $ immunedb_metadata --use-filenames


This creates a ``metadata.tsv`` file in ``/share/sequences`` in Docker or
``$HOME/immunedb_share/sequences`` on the host.

The ``--use-filenames`` flag is optional, and simply populates the
``sample_name`` field with the file names stripped of their ``.fasta`` or
``.fastq`` extension.

Editing the Metadata Sheet
--------------------------
On the host open the ``$HOME/immunedb_share/sequences`` file in Excel or your
favorite spreadsheet editor.  The headers included in the file are
**required**.  You may add additional headers as necessary for your dataset
(e.g. ``tissue``, ``cell_subset``, ``timepoint``) so long as they follow the
following rules:

* The headers must all be unique
* Each header may only contain *lowercase* letters, numbers, and underscores
* Each header must begin with a (lowercase) character
* Each header must not exceed 32 characters in length
* The *values* within each column cannot exceed 64 characters in length

.. note::

   When data is missing or not necessary in a field, leave it blank or set to
   ``NA``, ``N/A``, ``NULL``, or ``None`` (case-insensitive).

Running the Pipeline
====================
Much of the rest of the pipeline follows from the example pipeline's
:ref:`instance creation step <instance_creation>`.  To start, create a
database.  Here we'll call it ``my_db`` but you'll probably want to give it a
more descriptive name:


.. code-block:: bash
   :caption: Run in Docker

   $ immunedb_admin create my_db /share/configs

Then we'll identify the sequences.  For this process the germline genes must be
specified.  The germlines provided as FASTA files in the Docker image are:

* ``imgt_human_ighv`` & ``imgt_human_ighj``: Human B-cell heavy chains
* ``imgt_human_trav`` & ``imgt_human_traj``: Human T-cell α chains
* ``imgt_human_trbv`` & ``imgt_human_trbj``: Human T-cell β chains
* ``imgt_mouse_ighv`` & ``imgt_mouse_ighj``: Mouse B-cell heavy chains

For this segment we'll assume human B-cell heavy chains, but the process is the
same for any dataset:

.. code-block:: bash
   :caption: Run in Docker

   $ immunedb_identify /share/configs/my_db.json \
         /root/germlines/imgt_human_ighv.fasta \
         /root/germlines/imgt_human_ighj.fasta \
         /share/sequences
   $ immunedb_collapse /share/configs/my_db.json

Then we assign clones.  For B-cells we recommend:

.. code-block:: bash
   :caption: Run in Docker

   $ immunedb_clones /share/configs/my_db.json similarity

For T-cells we recommend:

.. code-block:: bash
   :caption: Run in Docker

   $ immunedb_clones /share/configs/my_db.json similarity --level nt \
         --min-similarity 1

If you have a mixed dataset, you can assign clones in different ways, filtering
on V-gene type.  For example:

.. code-block:: bash
   :caption: Run in Docker

   $ immunedb_clones /share/configs/my_db.json similarity --gene IGHV
   $ immunedb_clones /share/configs/my_db.json similarity --gene TCRB \
         --level nt --min-similarity 1


The last required step is to generate aggregate statistics:

.. code-block:: bash
   :caption: Run in Docker

    $ immunedb_clone_stats /share/configs/my_db.json
    $ immunedb_sample_stats /share/configs/my_db.json

For B-cells, you might want to generate lineages too.  The following excludes
mutations that only occur once.  ``immunedb_clone_trees`` has many other
parameters for filtering which you can view with the ``--help`` flag:

.. code-block:: bash
   :caption: Run in Docker

    $  immunedb_clone_trees /share/configs/my_db.json --min-mut-copies 2

Selection pressure can be run with the following.  This process is quite
time-consuming, even for small datasets:

.. code-block:: bash
   :caption: Run in Docker

    $ immunedb_clone_pressure /share/configs/my_db.json \
         /apps/baseline/Baseline_Main.r

Finally, the data should be available at http://localhost:8080/frontend/my_db.
