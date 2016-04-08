.. _docker_install:

Getting Started (Docker Compose)
================================

Dependency Installation
-----------------------
Using AIRRDB with docker requires that both Docker `Docker <http://docker.com>`_
and `Docker Compose <https://www.docker.com/products/docker-compose>`_ be
installed.  Please consult these websites to install them on your target system.

AIRRDB Docker Files
----------------------
Once docker is installed, run the following:

.. code-block:: bash

    $ mkdir ~/airrdb-docker
    $ cd ~/airrdb-docker
    $ wget https://github.com/arosenfeld/airrdb/blob/master/docker-compose.yml
    $ wget https://github.com/arosenfeld/airrdb/blob/master/run-docker.sh
    $ chmod +x run-docker.sh

This will download the Docker Compose file and a script to make running AIRRDB
easier.

Runnning AIRRDB
-------------
There are many ways to change the configuration of your Docker Compose AIRRDB
instance but to get started, within ``~/airrdb-docker`` you can run:

.. code-block:: bash

    $ ./run-docker.sh

This will do the following:

- Creates a directory ``~/airrdb/db`` which persistently stores the associated AIRRDB
  database.
- Creates a directory ``~/airrdb/data`` where you can place data in which you'd
  like to access from within the Docker instance.
- Starts a MariaDB (a MySQL drop-in replacement) database
- Starts REST service serving data from the database to a local port.
- Starts a ``airrdb-frontend`` process which serves a webpage in which to view
  data.
