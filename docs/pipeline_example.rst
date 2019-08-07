.. _pipeline_example:

Running the Example Pipeline
****************************
This page serves to familiarize new users with the basic flow of running the
ImmuneDB pipeline.  Example input FASTQ files are provided which contain human
B-cell heavy chain sequences.

Commands are listed as either being run in either the Docker container or on
the host.

To begin, run the Docker container :ref:`as documented
<running-the-container>`:

.. code-block:: bash
   :caption: Run on Host

    $ docker run -v $HOME/immunedb_share:/share \
         -p 8080:8080 -it arosenfeld/immunedb:v0.29.1


Metadata Specification
======================
Before ImmuneDB can be run, metadata must be specified for each input file.
For this example, one has already been created for you.  To learn how to create
a metadata file for your own data, see :ref:`Creating a Metadata Sheet`.


.. _instance_creation:


ImmuneDB Instance Creation
==========================
Next, we create a database for the data with:

.. code-block:: bash
   :caption: Run in Docker

    $ immunedb_admin create example_db /share/configs

This creates a new database named ``example_db`` and stores its configuration
in ``/share/configs/example_db.json``.


Identifying or Importing Sequences
==================================
The first step of the pipeline is to annotate sequences and store the resulting
data in the newly created database.  To do so, the ``immunedb_identify`` is
used.  It requires that V and J germline sequences be specified in two separate
FASTA files.  The Docker image provides Human & Mouse IGH, TRA, and TRB
germlines in ``$HOME/germlines``.

For this example, there are two provided input files in ``/example`` along with
the requisite ``metadata.tsv`` file which you can view with:

.. code-block:: bash
   :caption: Run in Docker

   $ ls /example

Given this, run the ``immunedb_identify`` command:

.. code-block:: bash
   :caption: Run in Docker

    $ immunedb_identify /share/configs/example_db.json \
         /root/germlines/imgt_human_ighv.fasta \
         /root/germlines/imgt_human_ighj.fasta \
         /example


Sequence Collapsing
===================
ImmuneDB determines the uniqueness of a sequence both at the sample and subject
level.  For the latter, ``immunedb_collapse`` is used to find sequences that are the
same except at positions that have an ``N``.  Thus, the sequences ``ATNN`` and
``ANCN`` would be collapsed.

To collapse sequences, run:

.. code-block:: bash
   :caption: Run in Docker

    $ immunedb_collapse /share/configs/example_db.json

Clonal Assignment
=================
After sequences are assigned V and J genes, they can be clustered into clones
based on CDR3 Amino Acid similarity with the ``immunedb_clones`` command.  This
takes a number of arguments which should be read before use.

There are three ways to create clones: based on CDR3 AA similarity, T-cell
exact CDR3 NT identity, and a lineage based method.  For this example we'll use
the similarity based method with default parameters:

.. code-block:: bash
   :caption: Run in Docker

    $ immunedb_clones /share/configs/example_db.json similarity

This will create clones where all sequences in a clone will have the same
V-gene, J-gene, and (by default) 85% CDR3 AA identity.

.. _stats_generation:

Statistics Generation
=====================
Two sets of statistics can be calculated in ImmuneDB:

- **Clone Statistics:** For each clone and sample combination, how many unique
  and total sequences appear as well as the mutations from the germline.
- **Sample Statistics:** Distribution of sequence and clone features on a
  per-sample basis, including V and J usage, nucleotides matching the germline,
  copy number, V length, and CDR3 length.  It calculates all of these with and
  without outliers, and including and excluding partial reads.

These are calculated with the ``immunedb_clone_stats`` and ``immunedb_sample_stats``
commands and must be run in that order.

.. code-block:: bash
   :caption: Run in Docker

    $ immunedb_clone_stats /share/configs/example_db.json
    $ immunedb_sample_stats /share/configs/example_db.json


Selection Pressure (Optional)
=============================

.. warning::
   Selection pressure calculations are time-consuming, so you can skip this
   step if time is limited.

Selection pressure of clones can be calculated with `Baseline
<http://selection.med.yale.edu/baseline/Archive>`_.  To do so run:

.. code-block:: bash
   :caption: Run in Docker

    $ immunedb_clone_pressure /share/configs/example_db.json \
         /apps/baseline/Baseline_Main.r

Note, this process is relatively slow and may take some time to complete.

.. _tree_generation:

Clone Trees (Optional)
======================
Lineage trees for clones is generated with the ``immunedb_clone_trees``
command.  The only currently supported method is neighbor-joining as provided
by `Clearcut <http://bioinformatics.hungry.com/clearcut>`_.

Among others, the ``--min-mut-copies`` parameter allows for mutations to be
omitted if they have not occurred at least a specified number of times.  This
can be useful to correct for sequencing error.


.. code-block:: bash
   :caption: Run in Docker

    $ immunedb_clone_trees /share/configs/example_db.json --min-mut-copies 2

Web Interface
=============
ImmuneDB has a web interface to interact with a database instance.  The Docker
container automatically makes this available at
http://localhost:8080/frontend/example_db

When you create more databases, simply replace `example_db` with the proper
databse name.
