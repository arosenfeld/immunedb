Data Analysis Pipeline
======================
The primary component of SLDB is its clonal identification pipeline which has
the capability to take as input raw sequences, determine likely V and J genes,
and finally group similar sequences into clones.

The pipeline is comprised of a number of steps which allows any portion of this
process to be replaced by another system.  For example, HighV-Quest could be
used for V and J assignment portion.  Further, the SLDB API allows developers to
integrate other tools into each step of the pipeline.

This page explains the basic workflow and assumes MySQL and SLDB are already
installed on the system.  It does not attempt to detail all the possible options
at each stage of the pipeline and users are encouraged to review the usage
documentation of each command.

SLDB Instance Creation
----------------------
It is assumed that the root user's username and password for MySQL is known.
To create a new SLDB instance, one can use ``sldb_admin``:

.. code-block:: bash

    $ sldb_admin create root DB_NAME CONFIG_DIR

Replacing ``DB_NAME`` with an appropriate database name and ``CONFIG_DIR`` with
a directory in which the database configuration will be stored will initialize
the instance.

.. note::

    You may use a user other than ``root`` so long as it has permissions to
    create databases, create users, and grant users permission to manipulate
    database in any way.

After running this, a database with the specified name will be created.  Further
a configuration file with the same name and a ``.json`` extension will be placed
in ``CONFIG_DIR``.  This configuration file will be the method of referencing
the database for the rest of the pipeline steps.

Data Preparation
----------------
Before running the SLDB pipeline, the input sequence data must be properly
structured.  All the FASTA files for analysis must be placed in a single
directory.  For example:

.. code-block:: bash

    $ ls sequence-data
    subjectABC_spleen.fasta
    subjectDEF_blood.fasta
    subjectXYZ_liver.fasta

SLDB needs some metadata about each of the FASTA files to process it.
Specifically, it **requires** the following information

- ``study_name``: The name of study for the sample (e.g. Lupus)
- ``sample_name``: The name of the sample.
- ``date``: The date the sample was acquired in the format YYYY-MM-DD.
- ``subject``: A unique identifier for the subject.  This must be unique to the
  entire SLDB instance as they are not contextual to the study.  Therefore if
  two studies use the same identifier for different subjects, they must be
  given new distinct identifiers.
- ``paired``: A ``true`` or ``false`` value (**not a string***) indicating if
  the reads in the sample are paired-end.

The following are optional for each file:

- ``subset``: The subset of the sample (e.g. Plasmablast, CD19+).  If none is
  specified, the field will be left blank.
- ``ig_class``: The isotype of the sample (e.g. IgE).
- ``tissue``: The tissue of the sample (e.g. Lung, Spleen).  If none is
  specified, the field will be left blank.
- ``disease``: Any disease(s) present in the subject when the sample was taken
  (e.g. Lupus).  If none is specified, the field will be left blank.
- ``lab``: The name of the lab sequencing the sample. If none is specified, the
  field will be left blank.
- ``experimenter``: The individual who sequenced the sample. If none is
  specified, the field will be left blank.

This information is specified in a ``metadata.json`` file which must be placed
in the same directory as the FASTA files.  The following is an example of such a
metadata file:

.. code-block:: json

    {
        "all": {
            "study_name": "Lupus",
            "paired": true
        },
        "subjectABC_spleen.fasta": {
            "subject": "ABC",
            "tissue": "Spleen",
            "date": "2015-09-13"
        },
        "subjectDEF_blood.fasta": {
            "subject": "DEF",
            "tissue": "Blood",
            "date": "2015-09-14"
        },
        "subjectXYZ_liver.fasta": {
            "subject": "XYZ",
            "tissue": "Liver",
            "date": "2015-09-15"
        }
    }


The ``all`` block applies the specified keys to all files in the directory (even
if they are not included in the metadata file).  If a key is specified both in
the ``all`` block and the block for a file, the value specified for the file is
used.

.. warning::
    It's advisable to not use terms like "None", "N/A", or an empty string to
    specify missing metadata.  Various portions of SLDB group information based
    on metadata, and will consider strings like these distinct from null
    metadata.

After creating the metadata file, the directory should look like:

.. code-block:: bash

    $ ls sequence-data
    metadata.json
    subjectABC_spleen.fasta
    subjectDEF_blood.fasta
    subjectXYZ_liver.fasta

Germline Files
--------------
SLDB requires that V and J germlines be specified in two separate FASTA files.
There are a number of restrictions on their format.

For V Germlines
^^^^^^^^^^^^^^^

- Genes must be in the format IGHVX*Y where X is the gene name and Y is the
  allele.  For example, IGHV1-18*01, IGHV5-a*03, and IGHV7-4-1*05 are all valid.
  However, IGHV1-18 and V1-18*01 are not.
- Germlines starting with gaps are excluded from alignment.
- Germlines must be IMGT gapped.
- V germlines must have have one of the following anchors with their last ``C``
  being the first base in the CDR3: ``D...Y[YCH]C``, ``Y[YHC]C`` or ``D.....C``.

For J Germlines
^^^^^^^^^^^^^^^
- Gene names follow the same rules as for V genes except they must start with
  ``IGHJ`` instead of ``IGHV``.
- There must be a fixed number of bases upstream of the CDR3 in all genes.

Sequence Identification
-----------------------
The first step of the pipeline is sequence identification.  Primarily this
assigns each sequence a V and J gene, but it also calculates statistics such as
how well the sequence matches the germline, if there is a probable insertion or
deletion, and how far into the CDR3 the V and J likely extend.

For identification a  FASTA file with IMGT aligned V germlines is required.
This can be downloaded from `IMGT's Gene-DB <http://imgt.org/genedb>`_ directly.

.. code-block:: bash

    $ sldb_identify /path/to/config.json /path/to/sequence-data-directory /path/to/v_germlines /path/to/j_germlines \
                    J_NTS_UPSTREAM_OF_CDR3 J_ANCHOR_SIZE J_MIN_ANCHOR_LEN

Where ``J_NTS_UPSTREAM_OF_CDR3`` are the fixed number of nucleotides in each
germline J gene upstream of the CDR3, J_ANCHOR_SIZE is the number of nucleotides
to use as an anchor, and J_MIN_ANCHOR_LEN dictates how many bases must match.
Their values for humans are 31, 18, and 12 respectively.  Graphically:

.. code-block:: bash

                                                |---- J_MIN_ANCHOR_LEN ----|
                                           |-------- J_ANCHOR_SIZE --------|
                 ...-- V --|-- CDR3 --|------ J_NTS_UPSTREAM_OF_CDR3 ------|
    j_germline:                 ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
    seq:         ...ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG


Local Alignment of Indel Sequences (Optional)
---------------------------------------------
After identification, certain sequences will be marked as being probable indels
(or misalignments).  To fix these, ``sldb_local_align`` can **optionally** be
used to properly gap sequences or germlines.  This process is inherently slow
and therefor may not be necessary in many cases.

.. code-block:: bash

    $ sldb_local_align /path/to/config.json /path/to/j_germlines \
                       J_NTS_UPSTREAM_OF_CDR3 J_ANCHOR_SIZE J_MIN_ANCHOR_LEN


Sequence Collapsing
------------------------------------
SLDB collapses sequences at two levels: the sample and the subject.  Collapsing
two sequences at a given level means that they are exactly the same when
excluding the positions where either sequence has an unknown base (``N``).
Thus, the sequences ``ATNN`` and ``ANCN`` would be collapsed.

This process is has been written in C rather than Python due to its
computational complexity.  This fact is transparent to the user, however.

To collapse sequences, run:

.. code-block:: bash

    $ sldb_collapse /path/to/config.json

The optional ``--subject-ids`` flag can specify that only samples from certain
subjects should be collapsed.

Clonal Assignment
-----------------
After sequences are assigned V and J genes, they can be clustered into clones
based on CDR3 Amino Acid similarity with the ``sldb_clones`` command.  This
takes a number of arguments which should be read before use.

A basic example of clonal assignment, not using all possible arguments:

.. code-block:: bash

    $ sldb_clones /path/to/config.json

.. _stats_generation:

Statistics Generation
---------------------
Two sets of statistics can be calculated in SLDB:

- **Sample Statistics:** Distribution of sequence and clone features on a
  per-sample basis, including V and J usage, nucleotides matching the germline,
  copy number, V length, and CDR3 length.  It calculates all of these with and
  without outliers, and including and excluding partial reads.
- **Clone Statistics:** For each clone and sample combination, how many unique
  and total sequences appear, mutations from the germline, and selection
  pressure.

These are calculated with the ``sldb_sample_stats`` and ``sldb_clone_stats``
commands.

For sample statistics there are only a few optional arguments which should be
reviewed.  In general, however, the command is issued to calculate statistics
for samples which do not already have them:

.. code-block:: bash

    $ sldb_sample_stats /path/to/config.json

Clone statistics require the path to the `Baseline
<http://selection.med.yale.edu/baseline/Archive>`_ main script.

.. code-block:: bash

    $ sldb_clone_stats /path/to/config.json /path/to/Baseline_Main.r

.. _tree_generation:

Clone Trees
-----------
Lineage trees for clones is generated with the ``sldb_clone_trees`` command.  The
only currently supported method is neighbor-joining as provided by `Clearcut
<http://bioinformatics.hungry.com/clearcut>`_.  Among others, the ``min-count``
parameter allows for mutations to be omitted if they have not occurred at least
a specified number of times.  This can be useful to correct for sequencing
error.


.. code-block:: bash

    $ sldb_clone_trees /path/to/config.json /path/to/clearcut --min-count 2

.. _supplemental_tools:

Supplemental Tools
------------------
In addition to the aforementioned pipeline commands, SLDB provides a number of
other commands.

sldb_hvquest
^^^^^^^^^^^^
This command can be used in place of ``sldb_identify`` to assign sequences V and
J genes from `HighV-Quest <http://www.imgt.org/HighV-QUEST>`_ output.  Since
there is no metadata file, all fields (e.g. subject, date, tissue) must be
manually specified.

Importing requires only two of the files output by HighV-Quest: the summary and
gapped nucleotides.

An example call to this command with only the required metadata:

.. code-block:: bash

    $ sldb_hvquest /path/to/config.json /path/to/summary_file \
        /path/to/gapped_nt_file /path/to/v_germlines STUDY_NAME SAMPLE_NAME
        READ_TYPE SUBJECT DATE

.. warning::
    SLDB may not be able to process some sequences from HighV-Quest, especially
    if it assigned a null CDR3.  Further, if the ``--v-ties`` flag is specified
    and the tied germline cannot be properly aligned to a sequence, it will be
    considered a no-result.


sldb_rest
^^^^^^^^^
SLDB has a RESTful API that allows for language agnostic querying.  This is
provided by the ``sldb_rest`` command.  It is specifically designed to provide
the required calls for the associated `web-app
<https://github.com/arosenfeld/simlab-web-database>`_.

It requires Haskell and the `diversity package
<https://hackage.haskell.org/package/diversity>`_.

To run on port 3000:

.. code-block::

    $ sldb_rest /path/to/config.json /path/to/diversity -p 3000
