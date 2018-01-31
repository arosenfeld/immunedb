Command Line Reference
**********************

This page provides brief descriptions of all the ImmuneDB commands and the
default help text for each.  For a step-by-step guide to running the pipeline,
refer to :doc:`pipeline`.

immunedb_admin
==============
Manages ImmuneDB instances.

.. program-output:: bin/immunedb_admin --help

immunedb_admin create
---------------------
Creates a new ImmuneDB instance.

.. program-output:: bin/immunedb_admin create --help

immunedb_admin delete
---------------------
Deletes an ImmuneDB instance.

.. program-output:: bin/immunedb_admin delete --help

immunedb_admin backup
---------------------
Performs a SQL dump from a ImmuneDB instance, backing up the contents.

.. program-output:: bin/immunedb_admin backup --help

immunedb_admin restore
----------------------
Restores a SQL dump from ``immunedb_admin backup``.

.. program-output:: bin/immunedb_admin restore --help

....

immunedb_clone_import
=====================
Imports custom clonal assignments from an exported sequence template.

To use, run ``immunedb_clone_import ... --action export``, fill in the
``clone_id`` column for each sequence, then run ``immunedb_clone_import ...
--action import``.  All sequences assigned to a given clone must have the same
V-gene, J-gene, and CDR3 length, otherwise an error will be raised.

.. program-output:: bin/immunedb_clone_import --help

....

immunedb_clone_pressure
=======================

.. note::

    For selection pressure calculations, `BASELINe
    <http://selection.med.yale.edu/baseline/>`_ must be installed.

Calculates the selection pressure acting on clonal sequences.

.. program-output:: bin/immunedb_clone_pressure --help

....

immunedb_clones
===============

Assigns sequences to clones using one of three methods.

.. program-output:: bin/immunedb_clones --help

....

immunedb_clones similarity
--------------------------

Groups sequences with the same subject, V-gene, J-gene, CDR3 length, and (by
default) 85% CDR3 AA similarity into clones.

.. program-output:: bin/immunedb_clones x similarity --help

immunedb_clones tcells
----------------------

Groups sequences with the same subject, V-gene, J-gene, and CDR3 NT sequence
into clones.

.. program-output:: bin/immunedb_clones x tcells --help

immunedb_clones lineage
-----------------------

.. note::

    To assign clones via the lineage tree method, `Clearcut
    <http://bioinformatics.hungry.com/clearcut/>`_ must be installed.

Creates lineages out of all sequences with the same subject, V-gene, J-gene,
and CDR3 length.  Then, the lineage is split along branches where the aggregate
number of mutations is at least ``--mut-cuttoff`` (default 4).

.. program-output:: bin/immunedb_clones x lineage --help

....

immunedb_clone_stats
====================
Aggregates statistics about clones for quicker, easier bulk querying.

.. program-output:: bin/immunedb_clone_stats --help

....

immunedb_clone_trees
====================
.. note::

    To create lineage trees, `Clearcut
    <http://bioinformatics.hungry.com/clearcut/>`_ must be installed.

Creates a lineage tree for each clone using Neighbor Joining.

.. program-output:: bin/immunedb_clone_trees --help

....

immunedb_collapse
=================
Collapses identical sequences across all samples in each subject.

.. program-output:: bin/immunedb_collapse --help

....

immunedb_export
===============
Exports data from ImmuneDB into various formats

.. program-output:: bin/immunedb_export --help

immunedb_export changeo
-----------------------

.. program-output:: bin/immunedb_export x changeo --help


immunedb_export genbank
-----------------------

.. program-output:: bin/immunedb_export x genbank --help


immunedb_export vdjtools
-----------------------

.. program-output:: bin/immunedb_export x vdjtools --help

....

immunedb_genotype
=================

.. note::

    To genotype subjects, `TIgGER <https://tigger.readthedocs.io>`_ must be
    installed.

Runs genotyping on a database that was generated with ``immunedb_identify ...
--genotyping``.

.. program-output:: bin/immunedb_genotype --help

....

immunedb_identify
=================

Identifies V- and J-genes of sequences in FASTA/FASTQ files `using an `anchor method
<https://www.ncbi.nlm.nih.gov/pubmed/26529062>`_.

.. program-output:: bin/immunedb_identify --help

....

immunedb_import
===============

.. note::

    Importing from a delimited file is still considered in beta and may not
    work as intended.  Please report any bugs on github.

Imports sequence alignments and gene calls from a tab delimited file (by
default in IMGT format)

.. program-output:: bin/immunedb_import --help

....

immunedb_local_align
====================

.. note::

    To locally align sequences `Bowtie2
    <http://bowtie-bio.sourceforge.net/bowtie2>`_ must be installed.

Corrects sequences that were flagged as potential indels or unidentifiable by
``immunedb_identify``.  This can be a slow process for large datasets.

.. program-output:: bin/immunedb_local_align --help

....

immunedb_metadata
=================

Generates metadata for a set of FASTA/FASTQ files to use for ``immunedb_identify``.

.. program-output:: bin/immunedb_metadata --help

....

immunedb_rest
=============

Starts a REST API server for ImmuneDB.  This can be used for any purpose, but
is designed to provide data for the `ImmuneDB Frontend
<github.com/arosenfeld/immunedb-frontend>`_.

.. program-output:: bin/immunedb_rest --help

....

immunedb_sample_stats
=====================

Calculates aggregate statistics for samples for faster querying.

.. program-output:: bin/immunedb_sample_stats --help

....

immunedb_sql
============

Starts an interactive MySQL session for a given ImmuneDB instance.  This is
simply a wrapper around the ``mysql`` command that passes information from a
configuration file.

.. program-output:: bin/immunedb_sql --help
