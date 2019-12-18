ImmuneDB |travisci| |docs| |codecov| |pypi| |docker|
====================================================
.. |travisci| image:: https://img.shields.io/travis/arosenfeld/immunedb/master.svg
  :target: https://travis-ci.org/arosenfeld/immunedb
.. |codecov| image:: https://img.shields.io/codecov/c/github/arosenfeld/immunedb.svg
  :target: https://codecov.io/gh/arosenfeld/immunedb
.. |docker| image:: https://img.shields.io/docker/pulls/arosenfeld/immunedb.svg
  :target: https://hub.docker.com/r/arosenfeld/immunedb
.. |pypi| image:: https://img.shields.io/pypi/v/immunedb.svg
  :target: https://pypi.python.org/pypi/ImmuneDB
.. |docs| image:: https://readthedocs.org/projects/immunedb/badge/?version=latest
   :target: https://immunedb.readthedocs.io/en/latest/?badge=latest

**ImmuneDB** is a database-backed system to analyze and store large amounts
(terabytes) of high-throughput B-cell receptor (BCR) and T-cell receptor (TCR)
data.  Although it can be used as a stand-alone package for comprehensive
repertoire profiling, ImmuneDB excels at acting as a central data store and
interface between other tools such as `IgBLAST
<https://ncbi.github.io/igblast>`_, the `Immcantation Framework
<http://immcantation.com>`_, `MiXCR <https://mixcr.readthedocs.io>`_, and
`VDJtools <https://vdjtools-doc.readthedocs.io>`_ via `AIRR compliant
<http://docs.airr-community.org/en/latest/resources/support.html#rearrangement-schema>`_
importing and exporting routines.

Feature Highlights
------------------

* **Relational storage of repertoire data**: Sequences, annotations, clones,
  lineages, and statistics are all stored in a relational database to
  promote consistent formatting and easy querying.
* **Consolidated metadata**: Custom study, experiment, and replicate metadata
  is stored alongside your sequencing data in a non-redundant format to avoid
  inconsistencies and errors over the life of your study.
* **Web interface**: ImmuneDB provides a built-in web interface for interactive
  exploration of data.
* **Interoperability**: With AIRR compliant input and output methods, ImmuneDB
  can interface with other software in the AIRR ecosystem.  Other output
  formats include Change-O and VDJtools.
* **Proven reliability**: ImmuneDB is used by multiple labs to manage terabytes
  of data comprised of billions of sequences in dozens of projects.


Quick Start
-----------
To get started immediately, please see the :doc:`Docker installation
instructions <install_docker>`.

.. toctree::
    :maxdepth: 2
    :hidden:
    :caption: Getting Started

    Introduction <self>

    install_docker
    install_local

.. toctree::
    :maxdepth: 2
    :hidden:
    :caption: Running the Pipeline

    pipeline_example
    pipeline_full
    modifying
    background

.. toctree::
    :maxdepth: 2
    :hidden:
    :caption: Using a Database

    exporting
    querying
    api

.. toctree::
    :maxdepth: 2
    :hidden:
    :caption: Other

    referencing
    related
