ImmuneDB |travisci| |codecov| |pypi| |docker|
=============================================
.. |travisci| image:: https://img.shields.io/travis/arosenfeld/immunedb/master.svg
  :target: https://travis-ci.org/arosenfeld/immunedb
.. |codecov| image:: https://img.shields.io/codecov/c/github/arosenfeld/immunedb.svg
  :target: https://codecov.io/gh/arosenfeld/immunedb
.. |docker| image:: https://img.shields.io/docker/pulls/arosenfeld/immunedb.svg
  :target: https://hub.docker.com/r/arosenfeld/immunedb
.. |pypi| image:: https://img.shields.io/pypi/v/immunedb.svg
  :target: https://pypi.python.org/pypi/ImmuneDB

**ImmuneDB** is a system to facilitate efficient storage and analysis
of high-throughput B- and T-cell sequence data.  It provides V- and J-gene
identification, clonal assignment, lineage construction, selection pressure
calculation, and thorough exporting functionality.

It also provides an intuitive and useful web interface a demo of which you can
see `here <http://immunedb.com/demo>`_.

Quick Start
-----------
To get started immediately, please see the :doc:`Docker installation
instructions <install_docker>`.

More Information
----------------

ImmuneDB is comprised of two GitHub repositories: Python analysis tools
(`arosenfeld/immunedb <https://github.com/arosenfeld/immunedb>`_) and a web
interface (`arosenfeld/immunedb-frontend
<https://github.com/arosenfeld/immunedb-frontend>`_)

The system aims to:

- **Reduce ad-hoc scripting:** Data analysis performed on an ad-hoc basis with
  custom scripts and data formats is error-prone and leads to inconsistencies.
  ImmuneDB provides a standardized analysis platform, performing many common tasks
  automatically.

- **Minimize flat-files:** Flat files are currently the standard method of
  data exchange in the biological sciences.  There are a myriad of drawbacks
  when using these including a lack of `referential integrity
  <http://en.wikipedia.org/wiki/Referential_integrity>`_, unclear `provenance
  <http://en.wikipedia.org/wiki/Provenance#Data_provenance>`_, and
  non-standardized formats.

  ImmuneDB attempts to reduce the need for flat files, through the use of an
  industry-leading database, MySQL.  When data must be exchanged as a flat-file,
  many export options, including FASTA and tab-delineation, are available.

- **Interoperate with existing tools:** ImmuneDB integrates tools from other
  researchers to provide features such as lineage construction, genotyping, and
  selection pressure calculations.  Further, ImmuneDB can import and export in
  a variety of common formats, making it compatible with the larger AIRR
  ecosystem of tools.

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

.. toctree::
    :maxdepth: 2
    :hidden:
    :caption: Querying an Existing Database

    api
    database

.. toctree::
    :maxdepth: 2
    :hidden:
    :caption: Other

    cli
    referencing
    related
