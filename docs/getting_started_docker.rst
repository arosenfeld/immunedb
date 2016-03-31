.. _docker_install:

Getting Started (Docker Compose)
================================

Dependency Installation
-----------------------
Using SLDB with docker requires that both Docker `Docker <http://docker.com>`_
and `Docker Compose <https://www.docker.com/products/docker-compose>`_ be
installed.  Please consult these websites to install them on your target system.

SLDB Docker Files
----------------------
Once docker is installed, run the following:

.. code-block:: bash

    $ mkdir ~/sldb-docker
    $ cd ~/sldb-docker
    $ wget https://github.com/arosenfeld/sldb/blob/master/docker-compose.yml
    $ wget https://github.com/arosenfeld/sldb/blob/master/run-docker.sh
    $ chmod +x run-docker.sh

This will download the Docker Compose file and a script to make running SLDB
easier.

Runnning SLDB
-------------
There are many ways to change the configuration of your Docker Compose SLDB
instance but to get started, within ``~/sldb-docker`` you can run:

.. code-block:: bash

    $ ./run-docker.sh

This will do the following:

- Creates a directory ``~/sldb/db`` which persistently stores the associated SLDB
  database.
- Creates a directory ``~/sldb/data`` where you can place data in which you'd
  like to access from within the Docker instance.
- Starts a MariaDB (a MySQL drop-in replacement) database
- Starts REST service serving data from the database to a local port.
- Starts a ``sldb-frontend`` process which serves a webpage in which to view
  data.
