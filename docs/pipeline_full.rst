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

Creating a Template Metadata Sheet
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
