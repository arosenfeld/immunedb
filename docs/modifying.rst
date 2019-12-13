.. _modifying:

Modifying the Database
======================
Databases can be modified in various ways using the ``immunedb_modify``
command.

Appending New Data
------------------
Adding new samples to a database is simply running the steps in
:ref:`pipeline_full` just on the new FASTA/FASTQ or AIRR files.  Effort has
been made to reduce the amount of information that needs to be recomputed when
samples are added.  However, after new samples are added all affected subjects
will be entirely re-collapsed and clones will be recalculated.


Changing Metadata
-----------------
Metadata specified when initially populating ImmuneDB via importing or
identification can be updated in two steps.  First, export the metadata
currently in the database with:

.. code-block:: bash

    $ immunedb_export PATH_TO_CONFIG samples --for-update

This will generate a ``samples.tsv`` file which can by modified.  Headers and
values can be changed, deleted, or added.

.. note::

    Note that changing the subject of any sample will require steps after and
    including ``immunedb_collapse`` to be re-run.

After modifying the metadata, update the database with:

.. code-block:: bash

    $ immunedb_modify PATH_TO_CONFIG update-metadata samples.tsv


Combining Samples
-----------------
.. warning::

    You cannot collapse samples from multiple subjects.  If that functionality
    is desired, first modify the metadata to set the subject for each sample to
    be the same with ``update-metadata``, and then run ``combine-samples``.

One assumption ImmuneDB makes is that each sample is a *biological replicate*
in that no one cell has its BCR/TCR sequence in more than one sample.  If you
have *technical replicates*, multiple independent sequencing runs of the same
same biological replicate, they should be combined into one ImmuneDB-sample
each.  To do so, add a metadata field to the database as described in
:ref:`Changing Metadata` where all technical replicates from the same
biological replicate have the same value.

For example, if we have the following samples where each sample has two
technical replicates:

================    =======
sample              subject
================    =======
biorep1_techrep1    S1
biorep1_techrep2    S1
biorep2_techrep1    S1
biorep2_techrep2    S1
================    =======

You would update the metadata to be:

================    =======     ========
sample              subject     collapse
================    =======     ========
biorep1_techrep1    S1          first_sample
biorep1_techrep2    S1          first_sample
biorep2_techrep1    S1          second_sample
biorep2_techrep2    S1          second_sample
================    =======     ========

And then run:

.. code-block:: bash

    $ immunedb_modify PATH_TO_CONFIG combine-samples collapse

This will result in the four replicates being collapsed into two, using the
``collapse`` field as the new name for each:

================    =======
sample              subject
================    =======
first_sample        S1
second_sample       S1
================    =======

Note the header ``collapse`` can have any value you want so long as it's passed
to ``immunedb_modify``.  Further, the values in that column can be arbitrary
but will be used as the new name of the samples after collapsing.

Deleting Samples
----------------
The following command can be used to delete samples by ID:


.. code-block:: bash

    $ immunedb_modify PATH_TO_CONFIG delete-samples [sample_ids]

Note that deleting samples will require the subject to be re-analyzed by
running all pipeline steps after and including ``immunedb_collapse``.
