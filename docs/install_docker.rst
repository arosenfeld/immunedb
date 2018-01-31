.. _docker_install:

Installing with Docker
======================

.. note::

    This installation method is still in development and may not work
    perfectly.  Please report any issues you encounter on github.

Dependency Installation
-----------------------
Using ImmuneDB with Docker requires that both Docker `Docker <http://docker.com>`_
and `Docker Compose <https://www.docker.com/products/docker-compose>`_ be
installed.  Please consult these websites to install them on your target system.

ImmuneDB Docker Files
---------------------
ImmuneDB comes with a script to simplify running it within Docker for simple
projects.  Once Docker and Docker Compose are installed, install the script and
related files with:

.. code-block:: bash

    $ wget -q -O - http://immunedb.com/docker-install | bash

Using the Helper Script
-----------------------
There are many ways to change the configuration of your Docker Compose ImmuneDB
instance but to get started, run:

.. code-block:: bash

    $ ./immunedb-docker.sh start

This will do the following:

- Creates a directory ``~/immunedb/db`` which persistently stores the associated
  ImmuneDB database.
- Creates a directory ``~/immunedb/data`` where you can place data in which you'd
  like to access from within the Docker instance.
- Starts a MariaDB (a MySQL drop-in replacement) database
- Starts REST service serving data from the database to a local port.
- Starts a ``immunedb-frontend`` process which serves a webpage in which to view
  data.

You can view the frontend website at ``http://localhost:8080``.  There will be
no data there since this is a fresh installation.  To launch a bash shell into
the ImmuneDB Docker image run:

.. code-block:: bash

    $ ./immunedb-docker.sh shell

To get data, such as germlines and sequences, into the container, copy them into
the host machine's ``~/immunedb/data`` directory.  It will then be available in the
Docker instance at ``~/data``.

**You're all set!** A configuration file for the database was automatically
created at ``~/configs/immunedb.json``. You can now follow the :ref:`pipeline
instructions <pipeline>`, skipping the database creation since that has been
done for you.

Place any input files you need in the created ``~/immunedb/data`` directory
which will then be accessible at ``/root/data`` in the container.

Advanced Configuration
----------------------
The ``immunedb-docker.sh`` script makes some assumptions about the host system:

- The REST API is running on the same machine at port 5000.
- That the URL used to access the frontend webpage is ``http://localhost:8080``.
- The database files should be stored at ``~/immunedb/db``.
- Shared data files are located on the host at ``~/immunedb/data``.

Any of these can be configured with the following environment variables:

- ``API_ENDPOINT``: Change the location from which the REST API should be
  accessed.  For example, if the hostname of the machine is ``sub.abc.com`` you
  can set this to be ``http://sub.abc.com:5000``.
- ``BASENAME``: Change the basename from which the website will be accessed.
  For example, if you're forwarding traffic from ``http://abc.com/example``,
  this would need to be set to ``example``.
- ``DB_VOLUME``: Change the location on the host where the database files should
  be stored.
- ``DATA_VOLUME``: Change the location on the host which is shared with the
  Docker image.  For example you may change it to something like
  ``~/sequence-data``.
- ``API_PORT``: Change the port on which the REST API is served.  This should
  match the port on ``API_ENDPOINT``.
- ``SERVE_PORT``: Change the port on which the website is served.

For example:

.. code-block:: bash

    $ API_ENDPOINT=http://sub.abc.com:5000 BASENAME=example SERVE_PORT=80 ./immunedb-docker.sh
