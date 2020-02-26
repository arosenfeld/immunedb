Running in the Background
*************************
After you have populated your ImmuneDB database(s), you may want to leave the
frontend web service running in the background.  To do so, you can start
ImmuneDB in detached mode with the following:

.. code-block:: bash

    $ docker run -v $HOME/immunedb_share:/share \
         -p 8080:8080 -e IMMUNEDB_DAEMON=1 -d=true \
         arosenfeld/immunedb:v0.29.10

If you want to stop the process in the future, get its process ID with

.. code-block:: bash

    $ docker ps


And then run:

.. code-block:: bash

    $ docker stop ID
