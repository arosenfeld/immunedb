Data Analysis Pipeline
======================
One of the primary components of SLDB is its clonal identification pipeline
which has the capability to take as input raw sequences, determine likely V and
J genes, and finally group similar sequences into clones.

The pipeline is comprised of a number of steps, however, which allows any
portion of this process to be replaced by another system.  For example,
HighV-Quest could be used for V and J assignment portion.  Further, the SLDB API
allows developers to integrate existing tools into each step of the pipeline.

This section will explain the basic workflow.  It assumes a set of FASTA files
in a single directory which are to be identified and assigned clones.

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
``HOST`` with the mySQL hostname (probably `localhost`), ``DATABASE`` with the
name of the master database, ``USER`` with the name of the user, and
``PASSWORD`` with the user's password.

.. code-block:: json

    {
        "host": "HOST",
        "database": "DATABASE",
        "username": "USER",
        "password": "PASSWORD"
    }

Then, create the file ``data.json`` with the same contents except replace
``DATABASE`` with the data database name.

Sequence Identification
-----------------------
The first step in the SLDB pipeline is sequence identification.  Primarily this
assigns each sequence a V and J gene, but it also calculates statistics such as
how well the sequence matches the germline, if there is a probable insertion or
deletion, and how far into the CDR3 the V and J likely extend.

First a metadata JSON file must be placed in the directory with the FASTA files.

.. code-block:: json

    {
        "all": {
            "key1": "value1",
            ...
            "keyn": "valuen"
        },
        "file1.fasta": {
            "key1": "value1",
            ...
            "keyn": "valuen"
        },
        "file2.fasta": {
            "key1": "value1",
            ...
            "keyn": "valuen"
        },
        ...
    }

For each file, the following keys **must** exist:

- ``read_type``: The type of read.  If the reads are paired, use ``R1+R2``,
  otherwise use ``R1`` or ``R2``.
- ``study_name``: The name of study for the sample (e.g. Lupus)
- ``subject``: A unique identifier for the subject.  If the SLDB instance has
  distinct subjects with the same identifier, they must be given different
  identifiers or they will be treated as the same subject.

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

The ``all`` block applies the specified keys to all files in the directory (even
if they are not included in the metadata file.  SLDB first looks in the block
for each file for each key.  If it is not found, only then does it check
``all``.

A small example is:

.. code-block:: json

    {
        "all": {
            "read_type": "R1+R2",
            "study_name": "Lupus"
        },
        "subjectABC_spleen.fasta": {
            "subject": "ABC",
            "tissue": "Spleen"
        },
        "subjectXYZ_liver.fasta": {
            "subject": "XYZ",
            "tissue": "Liver"
        }
    }

.. rubric:: Footnotes

.. [#clone_groups] With the exception of the ``clone_groups`` table which will
    potentially change.
