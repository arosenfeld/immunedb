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
Once docker is installed, clone the AIRRDB repository:

.. code-block:: bash

    $ git clone git@github.com:arosenfeld/airrdb.git ~/airrdb-git

In the resulting ``~/airrdb-git`` directory, there will be a ``run-docker.sh`` file
which will be used to launch the Docker instance.


Runnning AIRRDB
-------------
There are many ways to change the configuration of your Docker Compose AIRRDB
instance but to get started, within ``~/airrdb-git`` you can run:

.. code-block:: bash

    $ ./run-docker.sh

This will do the following:

- Creates a directory ``~/airrdb/db`` which persistently stores the associated
  AIRRDB database.
- Creates a directory ``~/airrdb/data`` where you can place data in which you'd
  like to access from within the Docker instance.
- Starts a MariaDB (a MySQL drop-in replacement) database
- Starts REST service serving data from the database to a local port.
- Starts a ``airrdb-frontend`` process which serves a webpage in which to view
  data.

You can view the frontend website at ``http://localhost:8080``.  There will be
no data there since this is a fresh installation.  To launch a bash shell into
the AIRRDB Docker image run the following

.. code-block:: bash

    $ docker ps
	CONTAINER ID        IMAGE                      COMMAND                  CREATED             STATUS              PORTS                    NAMES
	106fa8a740b3        arosenfeld/airrdb            "/bin/sh -c '/app/./w"   5 seconds ago       Up 4 seconds        0.0.0.0:5000->5000/tcp   airrdb_airrdb_1
	fe2895cf3be2        arosenfeld/airrdb-frontend   "npm run serve"          5 seconds ago       Up 4 seconds        0.0.0.0:8080->8080/tcp   airrdb_frontend_1
	18a58f85551e        mariadb                      "/docker-entrypoint.s"   6 seconds ago       Up 5 seconds        3306/tcp                 airrdb_mariadb_1

Note that there are three containers running; one for the front-end, one for the
database, and one for the AIRRDB Docker instance itself.  The first one (named
``airrdb_airrdb_1`` here, but may differ for each machine) is the AIRRDB Docker
instance.  Let's launch a bash shell into it:

.. code-block:: bash

    $ docker exec -i -t 106fa8a740b3 bash

Note you need not type the entire ``CONTAINER ID`` but just enough of it to
uniquely identify the container.  After running this, you will have a bash
shell into the instance.  AIRRDB is installed here.

To get data, such as germlines and sequences, into the container, copy them into
the host machine's ``~/airrdb/data`` directory.  It will then be available in the
Docker instance at ``~/data``.

Advanced Configuration
----------------------
The ``run-docker.sh`` script makes some assumptions about the host system:

- The REST API is running on the same machine at port 5000.
- That the URL used to access the frontend webpage is ``http://localhost:8080``.
- The database files should be stored at ``~/airrdb/db``.
- Shared data files are located on the host at ``~/airrdb/data``.

Any of these can be configured with the following environment variables:

- ``API_ENDPOINT``: Change the location from which the REST API should be
  accessed.  For example, if the hostname of the machine is ``sub.abc.com`` you
  can set this to be ``http://sub.abc.com:5000``.
- ``BASE_URL``: Change the location from which the website is accessed.  For
  example ``http://sub.abc.com``
- ``DB_VOLUME``: Change the location on the host where the database files should
  be stored.
- ``DATA_VOLUME``: Change the location on the host which is shared with the
  Docker image.  For example you may change it to something like
  ``~/sequence-data``.
- ``API_PORT``: Change the port on which the REST API is served.  This should
  match the port on ``API_ENDPOINT``.
- ``SERVE_PORT``: Change the port on which the website is served.  This should
  match the port on ``BASE_URL``.  If no port is specified in ``BASE_URL`` set
  this to 80 (default HTTP).

For example:

.. code-block:: bash

    $ API_ENDPOINT=http://sub.abc.com:5000 BASE_URL=http://sub.abc.com SERVE_PORT=80 ./run-docker.sh
