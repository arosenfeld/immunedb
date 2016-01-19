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
    the database in any way.

After running this, a database with the specified name will be created.  Further
a configuration file with the same name and a ``.json`` extension will be placed
in ``CONFIG_DIR``.  This configuration file will be the method of referencing
the database for the rest of the pipeline steps.

Data Preparation
----------------
Before running the SLDB pipeline, the input sequence data must be properly
structured.  Sequences must be separated into one file per sample.  That is,
sequences in the same file must be from the sequencing run or, conversely, that
sequences in different files could not have originated from the same cell.  This
is required for SLDB to properly count the number of unique sequences.

For example, a directory of FASTA files may look like this:

.. code-block:: bash

    $ ls sequence-data
    subjectABC_spleen.fasta
    subjectDEF_blood.fasta
    subjectXYZ_liver.fasta

SLDB needs some metadata about each of the FASTA files to process it.
Specifically, it **requires** the following information (the maximum number of
characters, when applicable, is shown in parenthesis):

- ``sample_name`` (128): The name of the sample.
- ``study_name`` (128): The name of study the sample belongs to.
- ``date``: The date the sample was acquired in YYYY-MM-DD format.
- ``subject`` (64): A unique identifier for the subject.  This must be unique to
  the entire SLDB instance as they are not contextual to the study.  Therefore
  if two studies use the same identifier for different subjects, they must be
  given new distinct identifiers.
- ``paired``: A ``true`` or ``false`` value (**not a string**) indicating if
  the reads in the sample are paired-end.

The following are **optional** for each file:

- ``subset`` (128): The subset of the sample (e.g. Plasmablast, CD19+).  If none is
  specified, the field will be left blank.
- ``tissue`` (16): The tissue of the sample (e.g. Lung, Spleen).  If none is
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
            "sample_name": "ABC_SPLEEN",
            "subject": "ABC",
            "tissue": "Spleen",
            "date": "2015-09-13"
        },
        "subjectDEF_blood.fasta": {
            "sample_name": "DEF_BLOOD",
            "subject": "DEF",
            "tissue": "Blood",
            "date": "2015-09-14"
        },
        "subjectXYZ_liver.fasta": {
            "sample_name": "XYZ_LIVER",
            "subject": "XYZ",
            "tissue": "Liver",
            "date": "2015-09-15"
        }
    }


The ``all`` block applies the specified keys to all files.  If a key is
specified both in the ``all`` block and the block for a file, the value
specified for the file is used.

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
There are a number of restrictions on their format.  Most common germlines can
be downloaded from `IMGT's Gene-DB <http://imgt.org/genedb>`_ directly.

For V Germlines
^^^^^^^^^^^^^^^

- Genes must be in the format IGHVX*Y where X is the gene name and Y is the
  allele.  For example, IGHV1-18*01, IGHV5-a*03, and IGHV7-4-1*05 are all valid.
  However, IGHV1-18 and V1-18*01 are not.
- Germlines starting with gaps are excluded from alignment.
- Germlines must be IMGT gapped.
- V germlines must have have one of the following amino-acid anchors with the
  trailing ``C`` being the first residue in the CDR3: ``D...Y[YCH]C``,
  ``Y[YHC]C`` or ``D.....C``.  The ``.`` character represents any amino acid,
  and ``[YHC]`` indicates any one of ``Y``, ``H``, or ``C``.

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

.. code-block:: bash

    $ sldb_identify /path/to/config.json /path/to/v_germlines.fasta /path/to/j_germlines.fasta \
                    J_NTS_UPSTREAM_OF_CDR3 J_ANCHOR_SIZE J_MIN_ANCHOR_LEN /path/to/sequence-data-directory

Where ``J_NTS_UPSTREAM_OF_CDR3`` are the fixed number of nucleotides in each
germline J gene upstream of the CDR3, ``J_ANCHOR_SIZE`` is the number of nucleotides
to use as an anchor, and ``J_MIN_ANCHOR_LEN`` dictates how many bases must match.
**Their values for IMGT human germlines are 31, 18, and 12 respectively**.  When
using other germlines, these values may be different.  The regions are shown
graphically below:

.. code-block:: bash

                                                |---- J_MIN_ANCHOR_LEN ----|
                                           |-------- J_ANCHOR_SIZE --------|
                 ...-- V --|-- CDR3 --|------ J_NTS_UPSTREAM_OF_CDR3 ------|
    j_germline:                 ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
    seq:         ...ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG


Local Alignment of Indel Sequences (Optional)
---------------------------------------------
.. warning::
    This step is still considered in beta.  Some corner cases may not be
    properly handled, and quality information from FASTQ files will not be
    included in aligned sequences.  Use this only if you can tolerate the
    possibility of errors or inconsistencies.

After identification, certain sequences will be marked as being probable indels
(or misalignments).  To fix these, ``sldb_local_align`` can **optionally** be
used to properly gap sequences or germlines.  This process is inherently slow
and therefor may not be necessary in many cases.  To use, the `seq-align
<https://github.com/noporpoise/seq-align>`_ package must be built and the path
to the resulting `needleman_wunsch` binary passed to SLDB.

.. code-block:: bash

    $ sldb_local_align /path/to/config.json /path/to/needleman_wunsch /path/to/j_germlines \
                       J_NTS_UPSTREAM_OF_CDR3


Sequence Collapsing
------------------------------------
SLDB determines the uniqueness of a sequence both at the sample and subject
level.  For the latter, ``sldb_collapse`` is used to find sequences that are the
same except at positions that have an ``N``.  Thus, the sequences ``ATNN`` and
``ANCN`` would be collapsed.

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

- **Clone Statistics:** For each clone and sample combination, how many unique
  and total sequences appear as well as the mutations from the germline.
- **Sample Statistics:** Distribution of sequence and clone features on a
  per-sample basis, including V and J usage, nucleotides matching the germline,
  copy number, V length, and CDR3 length.  It calculates all of these with and
  without outliers, and including and excluding partial reads.

These are calculated with the ``sldb_clone_stats`` and ``sldb_sample_stats``
commands and must be run in that order.

.. code-block:: bash

    $ sldb_sample_stats /path/to/config.json
    $ sldb_clone_stats /path/to/config.json


Selection Pressure (Optional)
-----------------------------
Selection pressure of clones can be calculated with `Baseline
<http://selection.med.yale.edu/baseline/Archive>`_.  After installing, run:

.. code-block:: bash

    $ sldb_clone_pressure /path/to/config.json /path/to/Baseline_Main.r

This process is relatively slow and may take some time to complete.

.. _tree_generation:

Clone Trees (Optional)
----------------------
Lineage trees for clones is generated with the ``sldb_clone_trees`` command.  The
only currently supported method is neighbor-joining as provided by `Clearcut
<http://bioinformatics.hungry.com/clearcut>`_.  Among others, the ``min-count``
parameter allows for mutations to be omitted if they have not occurred at least
a specified number of times.  This can be useful to correct for sequencing
error.


.. code-block:: bash

    $ sldb_clone_trees /path/to/config.json /path/to/clearcut --min-count 2

.. _supplemental_tools:


Web Service (Optional)
----------------------
SLDB has a RESTful API that allows for language agnostic querying.  This is
provided by the ``sldb_rest`` command.  It is specifically designed to provide
the required calls for the associated `web-app
<https://github.com/arosenfeld/sldb-frontend>`_.

To run on port 3000 for example:

.. code-block:: bash

    $ sldb_rest /path/to/config.json -p 3000
