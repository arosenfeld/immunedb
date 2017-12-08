.. _pipeline:

Data Analysis Pipeline
**********************
The primary component of ImmuneDB is its clonal identification pipeline which has
the capability to take as input raw sequences, determine likely V and J genes,
and finally group similar sequences into clones.

The pipeline is comprised of a number of steps which allows any portion of this
process to be replaced by another system.  For example, HighV-Quest could be
used for V and J assignment portion.  Further, the ImmuneDB API allows developers to
integrate other tools into each step of the pipeline.

This page explains the basic workflow and assumes MySQL and ImmuneDB are already
installed on the system.  It does not attempt to detail all the possible options
at each stage of the pipeline and users are encouraged to review the usage
documentation of each command.


Quick Start
===========
The details of all ImmuneDB commands are below, however, the following is the
basic set of commands to run ImmuneDB.  This assumes you have the V-germlines
in ``imgt_human_v.fasta`` and J-germlines in ``imgt_human_j.fasta`` within the
current working directory.  It also assumes there is a set of FASTA/FASTQ files
in the current working directory.


.. code-block:: bash

    # Navigate to directory with FASTA/FASTQ files
    $ cd path/to/my_sequences
    # Make a template metadata file
    $ immunedb_metadata
    # Edit the resulting metadata.tsv file
    # Create a database
    $ immunedb_admin create example_db ~/configs
    # Annotate sequences with V/J genes
    $ immunedb_identify ~/configs/example_db.json imgt_human_v.fasta imgt_human_j.fasta .
    # Collapse sequences across replicates
    $ immunedb_collapse ~/configs/example_db.json
    # Assemble clones
    $ immunedb_clones ~/configs/example_db.json similarity
    # Generate clone statistics
    $ immunedb_clone_stats ~/configs/example_db.json
    # Generate sample statistics
    $ immunedb_sample_stats ~/configs/example_db.json

Data Preparation
================
Before running the ImmuneDB pipeline, the input sequence data must be properly
structured.  Sequences must be separated into one file per sample.  That is,
sequences in the same file must be from the sequencing run or, conversely, that
sequences in different files could not have originated from the same cell.  This
is required for ImmuneDB to properly count the number of unique sequences.

For example, a directory of FASTA files may look like this:

.. code-block:: bash

    $ ls
    subjectD001_spleen.fasta
    subjectD002_blood.fasta
    subjectD002_liver.fasta

ImmuneDB needs some metadata about each of the FASTA files to process it.
Specifically, it **requires** the following information (the maximum number of
characters, when applicable, is shown in parenthesis):

- ``sample_name`` (128): The name of the sample.
- ``study_name`` (128): The name of study the sample belongs to.
- ``date``: The date the sample was acquired in YYYY-MM-DD format.
- ``subject`` (64): A unique identifier for the subject.  This must be unique to
  the entire ImmuneDB instance as they are not contextual to the study.  Therefore
  if two studies use the same identifier for different subjects, they must be
  given new distinct identifiers.

The following are **optional** for each file:

- ``subset`` (128): The subset of the sample (e.g. Plasmablast, CD19+).  If none is
  specified, the field will be left blank.
- ``tissue`` (32): The tissue of the sample (e.g. Lung, Spleen).  If none is
  specified, the field will be left blank.
- ``ig_class`` (8): The isotype of the sample (e.g. IgE).
- ``disease`` (32): Any disease(s) present in the subject when the sample was taken
  (e.g. Lupus).  If none is specified, the field will be left blank.
- ``lab`` (128): The name of the lab sequencing the sample. If none is specified, the
  field will be left blank.
- ``experimenter`` (128): The individual who sequenced the sample. If none is
  specified, the field will be left blank.
- ``v_primer`` (32): An arbitrary string indicating the V-gene primer used.
- ``j_primer`` (32): An arbitrary string indicating the J-gene primer used.

This information is specified in a ``metadata.tsv`` file which must be placed in
the same directory as the FASTA files.  A template for this file can be
generated based on FASTA files with:

A template for the files above could be generated with the following command,
while in the same directory:

.. code-block:: bash

    $ immunedb_metadata


The following is an example of such a metadata file with some information filled
in:

======================== ============ ============ ========== ======= ====== ======= ======= ======== ================= ======== ======== ========
file_name                study_name   sample_name  date       subject subset tissue  disease lab      experimenter      ig_class v_primer j_primer
======================== ============ ============ ========== ======= ====== ======= ======= ======== ================= ======== ======== ========
subjectD001_spleen.fasta B-cell Study D001_SPL     2015-09-13 D001    Naive  Spleen          Some lab Mr. Experimenter           Leader   J mix
subjectD002_blood.fasta  B-cell Study D002_BL      2015-09-14 D002    Naive  Blood           Some lab Mrs. Experimenter          Leader   J mix
subjectD002_liver.fasta  B-cell Study D002_Liver   2015-09-15 D003    Mature Liver           Some lab Mrs. Experimenter          FW1      J mix
======================== ============ ============ ========== ======= ====== ======= ======= ======== ================= ======== ======== ========

.. warning::
    Avoid using terms such as  "None", "N/A", or an empty string to specify
    missing metadata.  Various portions of ImmuneDB group information based on
    metadata, and will consider strings like these distinct from NULL metadata.

After creating the metadata file, the directory should look like:

.. code-block:: bash

    $ ls
    metadata.tsv
    subjectD001_spleen.fasta
    subjectD002_blood.fasta
    subjectD003_liver.fasta

Germline Files
--------------
ImmuneDB requires that V and J germlines be specified in two separate FASTA files.
There are a number of restrictions on their format.  Most common germlines can
be downloaded from `IMGT's Gene-DB <http://imgt.org/genedb>`_ directly.

V Germlines
^^^^^^^^^^^
- Genes must be in the format prefixX*Y or prefixX where X is the gene name and Y is the
  allele.  For example, IGHV1-18*01, TRBV5-a*03, and IGHV7-4-1 are all valid.
  However, IGHV4-34 is not.
- Germlines must be IMGT gapped.
- Germlines starting with gaps are excluded from alignment.
- For anchor identification,  ImmuneDB uses the V/J alignment method found in
  `PMID: 26529062`.  This requires V germlines to have have one of the
  following amino-acid anchors with the trailing ``C`` being the first residue
  in the CDR3: ``D...Y[YCH]C``, ``Y[YHC]C`` or ``D.....C``.  The ``.``
  character represents any amino acid, and ``[YHC]`` indicates any one of
  ``Y``, ``H``, or ``C``.  **Local alignment does not place these restrictions
  on germlines.**

J Germlines
^^^^^^^^^^^^^^^
- There must be a fixed number of bases upstream of the CDR3 in all genes.

Main Pipeline
=============
ImmuneDB Instance Creation
--------------------------
It is assumed that the root user's username and password for MySQL is known.
To create a new ImmuneDB instance, one can use ``immunedb_admin``:

.. code-block:: bash

    $ immunedb_admin create DB_NAME CONFIG_DIR

Replacing ``DB_NAME`` with an appropriate database name and ``CONFIG_DIR`` with
a directory in which the database configuration will be stored will initialize
the instance.

.. note::

    By default the root user is used to create the database.  You may use a user
    other than ``root`` with the ``--admin-user`` flag, so long as it has
    permissions to create databases, create users, and grant users permission to
    manipulate the database in any way.

After running this, a database with the specified name will be created.  Further
a configuration file with the same name and a ``.json`` extension will be placed
in ``CONFIG_DIR``.  This configuration file will be the method of referencing
the database for the rest of the pipeline steps.

Sequence Identification (Anchor method)
---------------------------------------
The first step of the pipeline is sequence identification.  Primarily this
assigns each sequence a V and J gene, but it also calculates statistics such as
how well the sequence matches the germline, if there is a probable insertion or
deletion, and how far into the CDR3 the V and J likely extend.

.. code-block:: bash

    $ immunedb_identify config.json v_germlines.fasta j_germlines.fasta \
        /path/to/sequence-data-directory

.. note::
    J-gene assignment requires three parameters, the number of nucleotides in
    the J after (upstream) of the CDR3, a conserved anchor size starting at the
    end of the J, and a minimum anchor length.  The J gene is searched for by
    using these anchors which are 31, 18 and 12 respectively in humans (and are
    the default values for ImmuneDB).  For other species, these values may need to
    be tweaked.  The regions are shown graphically below:

    .. code-block:: bash

                                               |---- J_MIN_ANCHOR_LEN ----|
                                               |-------- J_ANCHOR_SIZE --------|
                     ...-- V --|-- CDR3 --|------ J_NTS_UPSTREAM_OF_CDR3 ------|
        j_germline:                 ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
        seq:         ...ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG


Genotyping (Optional)
---------------------
.. warning::
    This step is still considered in beta.

ImmuneDB comes with a helper script to determine the genotype of subjects using
`TIgGER <https://tigger.readthedocs.io>`_.  This will determine which germline
V-genes are present in each subject, and if any contain novel mutations.  After
this determination, ImmuneDB can operate on the modified genotype FASTA file
for futher gene identification.

The basic process for this is to identify sequences at the allele level, export
sequences in Change-O format, run TIgGer to determine each subjects' genotype,
delete the originally identified sequences, and then re-run identification with
the new V-germlines.

.. code-block:: bash

    $ immunedb_admin create db_name ~/configs
    $ immunedb_identify ~/configs/db_name.json v_germlines.fasta j_germlines.fasta \
        /path/to/sequences --genotype
    $ immunedb_collapse ~/configs/db_name.json
    $ immunedb_export ~/configs/db_name.json changeo --min-subject-copies 2
    $ immunedb_genotype ~/configs/db_name.json v_germlines.fasta
    $ immunedb_admin delete ~/configs/db_name.json
    $ immunedb_admin create db_name ~/configs
    # For each subject
    $ immunedb_identify ~/configs/db_name.json SUBJECT.v_genotype.fasta j_germlines.fasta \
        /path/to/SUBJECT_sequence_data

Note in the final step (identifying sequences with the inferred genotype) you
must specify the sequences only associated with ``SUBJECT``.  This step must
then be repeated for each subject for which the genotype was inferred.


Local Alignment of Indel Sequences (Optional)
---------------------------------------------
.. warning::
    This step is still considered in beta.  Some corner cases may not be
    properly handled, and quality information from FASTQ files will not be
    included in aligned sequences.

After identification, certain sequences will be marked as being probable indels
(or misalignments).  To fix these, ``immunedb_local_align`` can **optionally** be
used to properly gap sequences or germlines.  It requires `bowtie2
<http://bowtie-bio.sourceforge.net/bowtie2>`_ to be installed and in your
``PATH``.

.. code-block:: bash

    $ immunedb_local_align config.json v_germlines.fasta j_germlines.fasta


Sequence Collapsing
------------------------------------
ImmuneDB determines the uniqueness of a sequence both at the sample and subject
level.  For the latter, ``immunedb_collapse`` is used to find sequences that are the
same except at positions that have an ``N``.  Thus, the sequences ``ATNN`` and
``ANCN`` would be collapsed.

This process is has been written in C rather than Python due to its
computational complexity.  This fact is transparent to the user, however.

To collapse sequences, run:

.. code-block:: bash

    $ immunedb_collapse config.json

The optional ``--subject-ids`` flag can specify that only samples from certain
subjects should be collapsed.

Clonal Assignment
-----------------
After sequences are assigned V and J genes, they can be clustered into clones
based on CDR3 Amino Acid similarity with the ``immunedb_clones`` command.  This
takes a number of arguments which should be read before use.

There are three ways to create clones: based on CDR3 AA similarity, T-cell
exact CDR3 NT identity, and a lineage based method.

Similarity Based
^^^^^^^^^^^^^^^^

A basic example of similarity-based clonal assignment, not using all possible
arguments:

.. code-block:: bash

    $ immunedb_clones config.json similarity

This will create clones where all sequences in a clone will have the same
V-gene, J-gene, and (by default) 85% CDR3 AA identity.

If you ran local-alignment on sequences, ImmuneDB can also associate clones
with insertions or deletions with a probable "parent" clone.  The parent clone
will have the same V-gene, J-gene, and CDR3 length.  Further, the CDR3 amino
acid sequences of the subclone will differ by no more than ``--min-similarity``
(default 85%).  This process can be enabled with ``--subclones``.

.. code-block:: bash

    $ immunedb_clones config.json --subclones similarity

T-cells
^^^^^^^

If your data is comprised of T-cell sequences, use the T-cell method:

.. code-block:: bash

    $ immunedb_clones config.json tcells

This will create clones from the sequences with the same V-gene, J-gene, and
identical CDR3 nucleotides.

Lineage Method
^^^^^^^^^^^^^^

.. warning::
    This clone assignment method is still considered in beta.

The lineage based method constructs a lineage for all sequences within
subjects that have the same V-gene, J-gene, and CDR3 NT length.  It then
splits the tree based on common mutations to create clones.

.. code-block:: bash

    $ immunedb_clones config.json lineage

Among other arguments, ``--mut-cuttoff`` (default 4) will determine how many
mutations must be in common for sequences to be placed in the same clone.

Importing Custom Assignments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you prefer to import your own clonal assignment, ImmuneDB allows you to
export sequences to a file which you can annotate with clone IDs.

.. code-block:: bash

    $ immunedb_clone_import config.json --action export sequences.tsv

This will generate a TSV file with all the unique sequences.  The last column,
``clone_id`` will be blank for all rows in the file.  To associate sequences
together as belonging to a clone, fill in the same value for each of their
``clone_id`` fields.  The value itself can be any string or integer, and only
serves as a unique identifier for each clone.

The sequences you assign to a given clone must belong to the same subject and
have the same V-gene, J-gene, and number of nucleotides in the CDR3.  Further,
changing any other values in the TSV file may lead to unpredictable results;
they are provided to give adequate information to external clonal assignment
programs.

Once the clones have been annotated:

.. code-block:: bash

    $ immunedb_clone_import config.json --action import sequences.tsv

.. _stats_generation:

Statistics Generation
---------------------
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

    $ immunedb_sample_stats config.json
    $ immunedb_clone_stats config.json


Selection Pressure (Optional)
-----------------------------
Selection pressure of clones can be calculated with `Baseline
<http://selection.med.yale.edu/baseline/Archive>`_.  After installing, run:

.. code-block:: bash

    $ immunedb_clone_pressure config.json /path/to/Baseline_Main.r

This process is relatively slow and may take some time to complete.

.. _tree_generation:

Clone Trees (Optional)
----------------------
Lineage trees for clones is generated with the ``immunedb_clone_trees`` command.  The
only currently supported method is neighbor-joining as provided by `Clearcut
<http://bioinformatics.hungry.com/clearcut>`_.  Among others, the ``min-count``
parameter allows for mutations to be omitted if they have not occurred at least
a specified number of times.  This can be useful to correct for sequencing
error.


.. code-block:: bash

    $ immunedb_clone_trees config.json /path/to/clearcut --min-count 2

.. _supplemental_tools:


Web Service (Optional)
----------------------
ImmuneDB has a RESTful API that allows for language agnostic querying.  This is
provided by the ``immunedb_rest`` command.  It is specifically designed to provide
the required calls for the associated `web-app
<https://github.com/arosenfeld/immunedb-frontend>`_.

To run on port 3000 for example:

.. code-block:: bash

    $ immunedb_rest config.json -p 3000

Optional Rollbar Support
^^^^^^^^^^^^^^^^^^^^^^^^
The server also has optional `Rollbar <https://rollbar.com/>`_ support, allowing
the database maintainer to monitor for errors.  Before using Rollbar you must
install its package with ``pip install rollbar`` and get a Rollbar token from
their website.  Then, you can use it with:

.. code-block:: bash

    $ immunedb_rest config.json --rollbar-token YOUR_TOKEN

There is also the optional ``--rollbar-env NAME`` parameter which allows you to
specify the environment name for Rollbar (defaults to ``develop``).
