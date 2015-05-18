Data Analysis Pipeline
======================
One of the primary components of SLDB is its clonal identification pipeline
which has the capability to take as input raw sequences, determine likely V and
J genes, and finally group similar sequences into clones.

The pipeline is comprised of a number of steps, however, which allows any
portion of this process to be replaced by another system.  For example,
HighV-Quest could be used for V and J assignment portion.  Further, the SLDB API
allows developers to integrate other tools into each step of the pipeline.

This page explains the basic workflow and assumes MySQL and SLDB are already
installed on the system.

SLDB Instance Creation
----------------------
To create a new SLDB instance, first two empty databases must be created, one
master and one data.  To do so run the following replacing ``USER`` with a
username that has privileges to create databases, ``MASTER`` with the name for
the master database, and ``DATA`` with the name for the data database.

.. code-block:: bash

    $ mysql -u USER -p -e "create database MASTER; create database DATA;"

Then two configuration files must be made for SLDB to access each of these
databases.  Create a file ``master.json`` with the following contents, replacing
``HOST`` with the MySQL hostname (probably `localhost`), ``DATABASE`` with the
name of the **master** database, ``USER`` with the name of the user, and
``PASSWORD`` with the user's password.

.. code-block:: json

    {
        "host": "HOST",
        "database": "DATABASE",
        "username": "USER",
        "password": "PASSWORD"
    }

Then, create the file ``data.json`` with the same contents except replace
``DATABASE`` with the **data** database name.

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

- ``read_type``: The type of read.  If the sequences are full reads including
  the entire V and J (as defined by IMGT germlines), use ``R1+R2``.  If the
  reads are partial from the 3' (J gene) end use ``R1``, and if the reads are
  from the 5' (V gene) end use ``R2``.
- ``study_name``: The name of study for the sample (e.g. Lupus)
- ``subject``: A unique identifier for the subject.  This must be unique to the
  entire SLDB instance as they are not contextual to the study.  Therefore if
  two studies use the same identifier for different subjects, they must be
  given new distinct identifiers.

The following are optional for each file:

- ``date``: The date the sample was acquired in the format YYYY-MM-DD.  If none
  is specified, the field is left blank.
- ``sample_name``: The name of the sample.  If none is specified, the basename
  (the filename prior to the first ".") will be used.
- ``subset``: The subset of the sample (e.g. Plasmablast, CD19+).  If none is
  specified, the field will be left blank.
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
            "read_type": "R1+R2",
            "study_name": "Lupus"
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
    Do not use terms like "None", "N/A", or an empty string to specify missing
    metadata.  Various portions of SLDB group information based on metadata, and
    will consider strings like these distinct from null metadata.

After creating the metadata file, the directory should look like:

.. code-block:: bash

    $ ls sequence-data
    metadata.json
    subjectABC_spleen.fasta
    subjectDEF_blood.fasta
    subjectXYZ_liver.fasta

Sequence Identification
-----------------------
The first step of the pipeline is sequence identification.  Primarily this
assigns each sequence a V and J gene, but it also calculates statistics such as
how well the sequence matches the germline, if there is a probable insertion or
deletion, and how far into the CDR3 the V and J likely extend.

For identification a  FASTA file with IMGT aligned V germlines is required.
This can be downloaded from `IMGT's Gene-DB <http://imgt.org/genedb>`_ directly.

To run identification, the ``sldb_identify`` command is used.  All SLDB commands
can be passed the ``--help`` flag to print the usage instructions.  The basic
usage for identification requires the master config, data config, the germline
FASTA file, and the path to the directory with the metadata and FASTA files to
identify:

.. code-block:: bash

    $ sldb_identify /path/to/master.json /path/to/data.json /path/to/sequence-data-directory /path/to/v_germlines

Sample & Subject Sequence Collapsing
------------------------------------
SLDB collapses sequences at three levels: the sample, subject, and clone.
Collapsing two sequences at a given level means that they are considered the
same.  This is not necessarily a simple task of checking for exactly matching
sequences because sequences may contain unknown bases, represented by the
character ``N``.

To do so, SLDB compares sequences pairwise, ignoring positions that either
sequences has an ``N`` with the ``sldb_collapse`` command.  This process is
highly optimized due to its computational complexity, and written in C rather
than Python.

Collapsing at the sample level must be done immediately after identification
with:

.. code-block:: bash

    $ sldb_collapse /path/to/master.json /path/to/data.json samples

The optional ``--sample-ids`` flag can specify that only certain samples should be
collapsed.

After this, the sequences must be collapsed to the subject level with a similar
command:

.. code-block:: bash

    $ sldb_collapse /path/to/master.json /path/to/data.json subjects

Similarly, the ``--subject-ids`` flag can be used to limit which subjects are
collapsed.

Clonal Assignment
-----------------
After sequences are assigned V and J genes, they can be clustered into clones
based on CDR3 Amino Acid similarity with the ``sldb_clones`` command.  This
takes a number of arguments which should be read before use.

A basic example of clonal assignment, not using all possible arguments:

.. code-block:: bash

    $ sldb_clones /path/to/master.json /path/to/data.json --similarity 65 -order

This will assign each sequence with at least 2 copies to a clone.  Additionally,
it will establish clone-groups in the master database which make associating
clones across versions simpler.


Clone Sequence Collapsing
-------------------------
After clones are created, sequences within them must further be collapsed with:

.. code-block:: bash

    $ sldb_collapse /path/to/master.json /path/to/data.json clones

As with subject collapsing ``--subject-ids`` can be used to limit the clones to
specified subjects.

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

    $ sldb_sample_stats /path/to/master.json /path/to/data.json

Clone statistics require the path to the `Baseline
<http://selection.med.yale.edu/baseline/Archive>`_ main script.

.. code-block:: bash

    $ sldb_clone_stats /path/to/master.json /path/to/data.json /path/to/Baseline_Main.r

.. _tree_generation:

Clone Trees
-----------
Lineage trees for clones is generated with the ``sldb_clone_tree`` command.  The
only currently supported method is neighbor-joining as provided by `Clearcut
<http://bioinformatics.hungry.com/clearcut>`_.  Among others, the ``min-count``
parameter allows for mutations to be omitted if they have not occurred at least
a specified number of times.  This can be useful to correct for sequencing
error.


.. code-block:: bash

    $ sldb_clone_tree /path/to/master.json /path/to/data.json /path/to/Baseline_Main.r nj /path/to/clearcut --min-count 2

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

    $ sldb_hvquest /path/to/master.json /path/to/data.json /path/to/summary_file \
        /path/to/gapped_nt_file /path/to/v_germlines STUDY_NAME SAMPLE_NAME
        READ_TYPE SUBJECT DATE

.. warning::
    SLDB may not be able to process some sequences from HighV-Quest, especially
    if it assigned a null CDR3.  Further, if the ``--v-ties`` flag is specified
    and the tied germline cannot be properly aligned to a sequence, it will be
    considered a no-result.

sldb_modify_clone
^^^^^^^^^^^^^^^^^
In some cases, manually changing clone attributes may be desired.  For example,
some individuals may have clones with insertions or deletions in their germline
which requires a change to the V gene.  This can be achieved with
``sldb_modify_clone``.

For example, to add three gaps at position 70 to clone 1234:

.. code-block:: bash

    $ sldb_modify_clone /path/to/master.json /path/to/data.json 1234 70,3 --v-name IGHV3-34*01_deletion

It is necessary to specify a new V gene name since SLDB assumes germlines with
the same name have the exact same sequence.

Other operations are possible with ``sldb_modify_clone`` which can be shown with
the ``--help`` flag.

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

    $ sldb_rest /path/to/master.json /path/to/data.json /path/to/diversity -p 3000
