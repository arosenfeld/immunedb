.. _pipeline_full:

Running the Pipeline on Your Data
*********************************
This page describes how to run the ImmuneDB pipeline on your own BCR/TCR data.
It is assumed that you've previously tried the :ref:`example pipeline
<pipeline_example>` and understand the basics of running commands in the Docker
container.

Like in the example, each code block has a header saying if the command should
be run on the host or in the Docker container.

Copying Your Sequence Data Into Docker
======================================
Unlike in the :ref:`example pipeline <pipeline_example>` where sequencing data
was provided, you'll need to copy your own FASTA/FASTQ sequencing data or
AIRR-formatted IgBLAST output into the Docker container.

To do so, on the host, we create a new directory in the shared directory
into which we'll copy your sequencing data.  Here we're calling it
``sequences`` but you'll probably want to choose a more descriptive name.
Replace ``PATH_TO_SEQUENCES`` with the path to your sequencing data.

.. code-block:: bash
   :caption: Run on Host

   $ mkdir -p $HOME/immunedb_share/input
   $ cp PATH_TO_SEQUENCES $HOME/immunedb_share/input


Running IgBLAST (optional)
==========================

.. note::

    If your data is already in AIRR-compliant IgBLAST format or you are
    planning on using the built in anchoring method, you can skip this step.

The following command will run IgBLAST on your files.  Valid values for species
and locus are:

* ``SPECIES``:``human``, ``mouse``
* ``LOCUS``: ``IGH``, ``IGL``, ``IGK``, ``TRA``, ``TRB``,

.. code-block:: bash

    $ run_igblast.sh SPECIES LOCUS /share/input /share/input

For consistency with the commands in the rest of this tutorial, we'll move the
new IgBLAST output files to ``/share/input`` and move the FASTA/FASTQ files to
``/share/sequences``.

.. code-block:: bash

    $ mkdir -p /share/sequences
    $ mv /share/input/*.fast[aq] /share/sequences

Creating a Metadata Sheet
=========================
Next, we'll use the ``immunedb_metadata`` command to create a template metadata
file for your sequencing data.  In the Docker container run:

.. code-block:: bash
   :caption: Run in Docker

   $ cd /share/input
   $ immunedb_metadata --use-filenames

.. note::

    This command expects the files to end in .fasta for FASTA, .fastq for
    FASTQ, or .tsv for AIRR.

This creates a ``metadata.tsv`` file in ``/share/input`` in Docker which is
linked to ``$HOME/immunedb_share/input`` on the host.

The ``--use-filenames`` flag is optional, and simply populates the
``sample_name`` field with the file names stripped of their extension.

Editing the Metadata Sheet
--------------------------
On the host open the metadata file in Excel or your favorite spreadsheet
editor.  The headers included in the file are **required**.  You may add
additional headers as necessary for your dataset (e.g. ``tissue``,
``cell_subset``, ``timepoint``) so long as they follow the following rules:

* The headers must all be unique
* Each header may only contain *lowercase* letters, numbers, and underscores
* Each header must begin with a (lowercase) character
* Each header must not exceed 32 characters in length
* The *values* within each column cannot exceed 64 characters in length

.. note::

   When data is missing or not necessary in a field, leave it blank or set to
   ``NA``, ``N/A``, ``NULL``, or ``None`` (case-insensitive).

Pipeline Steps
==============
Much of the rest of the pipeline follows from the example pipeline's
:ref:`instance creation step <instance_creation>`.  To start, create a
database.  Here we'll call it ``my_db`` but you'll probably want to give it a
more descriptive name:


.. code-block:: bash
   :caption: Run in Docker

   $ immunedb_admin create my_db /share/configs

Then we'll identify or import the sequences.  For this process the germline
genes must be specified.  The germlines are provided FASTA files in the Docker
image at ``/root/germlines``.

.. note::

    You can use your own germline files if you desire so long as they are IMGT
    gapped.

For this segment we'll assume human B-cell heavy chains, but the process is the
same for any dataset.  Depending on if you want to use IgBLAST input
(recommended) or the built-in annotation method the command will be one of the
following:

**Option 1: Importing from IgBLAST output (recommended)**:

.. code-block:: bash
    :caption: Run in Docker

    $ immunedb_import /share/configs/example_db.json airr \
         /root/germlines/igblast/human/IGHV.gapped.fasta \
         /root/germlines/igblast/human/IGHJ.gapped.fasta \
         /share/input

**Option 2: Using anchoring method**:

.. code-block:: bash
   :caption: Run in Docker

   $ immunedb_identify /share/configs/my_db.json \
         /root/germlines/anchor/human/IGHV.gapped.fasta \
         /root/germlines/anchor/human/IGHJ.gapped.fasta \
         /share/input

After importing or identifying sequences, continue running the pipeline from
here:

.. code-block:: bash
    :caption: Run in Docker

    $ immunedb_collapse /share/configs/my_db.json

Then we assign clones.  For B-cells we recommend:

.. code-block:: bash
   :caption: Run in Docker

   $ immunedb_clones /share/configs/my_db.json cluster

For T-cells we recommend:

.. code-block:: bash
   :caption: Run in Docker

   $ immunedb_clones /share/configs/my_db.json cluster --min-similarity 1

If you have a mixed dataset, you can assign clones in different ways, filtering
on V-gene type.  For example:

.. code-block:: bash
   :caption: Run in Docker

   $ immunedb_clones /share/configs/my_db.json cluster --gene IGHV
   $ immunedb_clones /share/configs/my_db.json cluster --gene TCRB \
         --min-similarity 1


The last required step is to generate aggregate statistics:

.. code-block:: bash
   :caption: Run in Docker

    $ immunedb_clone_stats /share/configs/my_db.json
    $ immunedb_sample_stats /share/configs/my_db.json

For B-cells, you might want to generate lineages too.  The following excludes
mutations that only occur once.  ``immunedb_clone_trees`` has many other
parameters for filtering which you can view with the ``--help`` flag or at
:ref:`clone_trees`.

.. code-block:: bash
   :caption: Run in Docker

    $  immunedb_clone_trees /share/configs/my_db.json --min-seq-copies 2

Selection pressure can be run with the following.  This process is quite
time-consuming, even for small datasets:

.. code-block:: bash
   :caption: Run in Docker

    $ immunedb_clone_pressure /share/configs/my_db.json \
         /apps/baseline/Baseline_Main.r

Finally, the data should be available at http://localhost:8080/frontend/my_db.

Analyzing Your Data
===================
After all the above steps are complete, you should have a fully populated
database, ready for analysis via :ref:`exporting`, :ref:`querying`, and the
:ref:`api`.
