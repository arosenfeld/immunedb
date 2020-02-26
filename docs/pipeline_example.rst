.. _pipeline_example:

Running the Example Pipeline
****************************
This page serves to familiarize new users with the basic process of running the
ImmuneDB pipeline.  Example input FASTQ files are provided which contain human
B-cell heavy chain sequences.

Commands are listed as either being run in either the Docker container or on
the host.  All ``immunedb_*`` commands have a ``--help`` flag which will show
all arguments and their descriptions.  It is recommended for each command you
run the help flag to see options not listed in this documentation.

To begin, run the Docker container :ref:`as documented
<running-the-container>`:

.. code-block:: bash
   :caption: Run on Host

    $ docker run -v $HOME/immunedb_share:/share \
         -p 8080:8080 -it arosenfeld/immunedb:v0.29.10


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
Data can be added to the new ImmuneDB database either by importing annotated
sequencing data in `AIRR format
<http://docs.airr-community.org/en/latest/datarep/rearrangements.html>`_, or
via a built-in gene assignment method based on `Zhang, et al., 2015
<https://www.ncbi.nlm.nih.gov/pubmed/26529062>`_.

For this example, there are two input FASTQ files in ``/example/fastq`` along
with an associated metadata file.  These will be used regardless of the method
you choose.  There are also germline files for human and mouse included.

Option 1: Importing from AIRR Files (Recommended)
-------------------------------------------------
First, IgBLAST needs to be run on the input files.  A small wrapper script is
provided in the Docker container.  It takes three parameters: a species, a
locus (in uppercase), an input directory with FASTA/FASTQ files, and an output
directory:

.. code-block:: bash

    $ run_igblast.sh human IGH /example/fastq /example/airr

To import these files, run the ``immunedb_import`` command:

.. code-block:: bash
    :caption: Run in Docker

    $ immunedb_import /share/configs/example_db.json airr \
         /root/germlines/igblast/human/IGHV.gapped.fasta \
         /root/germlines/igblast/human/IGHJ.gapped.fasta \
         /example/airr

Option 2: Annotating FASTA/FASTQ Files via Anchoring
----------------------------------------------------
Alternatively, if you'd prefer to use the built-in annotation method on
FASTA/FASTQ files, you can use the ``immunedb_identify`` command.  Note this
method is more sensitive to high mutation rates in the regions flanking the
CDR3.

.. code-block:: bash
   :caption: Run in Docker

    $ immunedb_identify /share/configs/example_db.json \
         /root/germlines/anchor/human/IGHV.gapped.fasta \
         /root/germlines/anchor/human/IGHJ.gapped.fasta \
         /example


Sequence Collapsing
===================
After data are imported or annotated on a sample-level basis, ImmuneDB
determines the subject-level unique sequences; that is, the set of unique
sequences across all samples in each subject.  Because sequences may contain
the ambiguous ``N`` symbol, the process is not trivial string equality
checking.  It is implemented in the ``immunedb_collapse`` command.

To collapse sequences, run:

.. code-block:: bash
   :caption: Run in Docker

    $ immunedb_collapse /share/configs/example_db.json

Clonal Assignment
=================
After collapsing unique sequences across each subject they can be grouped into
clones which are aggregations of sequences likely deriving from a common
progenitor cell.

ImmuneDB offers two clonal inference methods, *similarity* and *cluster*.  The
*cluster* method is recommended and documented here as it more robust than
*similarity*.

For both methods, clones are inferred in two steps: grouping sequences and then
merging similar clones.  Both steps are run together with the
``immunedb_clones`` command

By default, only sequences with a subject-level copy number greater
than 1 are included in clones.  This can be changed with the ``--min-copy``
parameter.

In the first step of clonal inference, sequences meeting the above copy number
criteria are hierarchically clustered together such that any two sequences in a
clone must (1) have the same CDR3 length and (2) share at least 85% amino-acid
similarity in the CDR3.  The similarity can be changed with ``--min-similarity
X`` parameter where X is the minimum similarity between 0 and 1.  If nucleotide
similarity should be used, ``--level nt`` can be passed.

.. note::

    For T-cells it is recommended the ``--min-similarity 1`` parameter be set
    but the ``--level`` parameter by left at the default amino-acid setting.
    Using both ``--min-similarity 1 --level nt`` may lead to the creation of
    spurious clones due to sequencing error.  Only pass both if you're quite
    certain your sequencing error has been eliminated (e.g. by barcoding).

After this step is complete, sequences have been assigned to clones.  In some
cases clones may have the same CDR3 nucleotide sequence but different gene
calls.  This can indeed occur biologically but frequently due to mutation and
sequencing error causing incorrect gene calls.

To rectify this, a second step in clonal inference is to collapse merge clones
that have the same CDR3 nucleotide sequences.  In cases where this occurs, the
highest copy clone absorbs the lower copy clones.  This second step can be
configured in two ways via the ``--reduce-difference`` flag.  Setting it to a
negative number (e.g. ``--reduce-difference -1``) disables the step entirely.
Setting it to a positive number (e.g. ``--reduce-difference 2``) will alter the
step's behavior to combine clones differing by at most that number of
nucleotides.  The default value is 0, so only clones with exactly the same CDR3
nucleotide sequences will be combined.


.. code-block:: bash
   :caption: Run in Docker

    $ immunedb_clones /share/configs/example_db.json cluster

.. _stats_generation:

Statistics Generation
=====================
Two sets of statistics can be calculated in ImmuneDB:

- **Clone Statistics:** For each clone and sample combination, statistics on
  the clone's size, mutation level, and top copy sequence
- **Sample Statistics:** Distribution of sequence and clone features on a
  per-sample basis, including gene usage, mutation level, copy number, CDR3
  length.

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

.. _clone_trees:

Clone Trees (Optional)
======================
Lineage trees for clones is generated with the ``immunedb_clone_trees``
command.  The only currently supported method is neighbor-joining as provided
by `Clearcut <http://bioinformatics.hungry.com/clearcut>`_.

There are many parameters that can be changed for tree construction:

* ``--min-seq-copies`` (default 0): The minimum number copy number required for
  a sequence to be included in the tree.
* ``--min-seq-samples`` (default 0): The minimum number samples in which a
  sequence must appear for it to be included in the tree.
* ``--min-mut-copies`` (default 0): The minimum number of copies in which a
  mutation must occur to be included in the tree.
* ``--min-mut-samples`` (default 0): The minimum number of samples in which a
  mutation must occur to be included in the tree.
* ``--exclude-stops`` (default ``False``): Exclude sequences with a stop codon.
* ``--full-seq`` (default ``False``): By default only the V-region of each
  sequence (the portion 5' of the CDR3) is included in the tree construction.
  Setting this flag will use the entire sequence.

Generally we recommend using ``--min-seq-copies 2``.

.. code-block:: bash
   :caption: Run in Docker

    $ immunedb_clone_trees /share/configs/example_db.json --min-seq-copies 2

Web Interface
=============
ImmuneDB has a web interface to interact with a database instance.  The Docker
container automatically makes this available at
http://localhost:8080/frontend/example_db

When you create more databases, simply replace `example_db` with the proper
database name.

Next Steps
==========
Now that the basic workflow has been covered, instructions to run ImmuneDB on
your own data can be found at :ref:`pipeline_full`.
