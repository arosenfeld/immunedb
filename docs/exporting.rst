.. _exporting:

Exporting Data to Files
=======================
You can use the ``immunedb_export`` command to export your data in a variety of
formats.

Exporting Samples
-----------------
To export samples statistics run the command:

.. code-block:: bash

    $ immunedb_export PATH_TO_CONFIG samples

After completion, a TSV file ``samples.tsv`` will be written with the following
headers, one line per sample:

================================= ===========
Field                             Description
================================= ===========
``id``                            Unique numeric sample identifier
``name``                          Name given to the sample
``subject``                       Subject from which the sample originated
``input_sequences``               Reads input into ImmuneDB
``identified``                    Reads successfully annotated
``in_frame``                      Reads in-frame
``stops``                         Reads with stop codons
``functional``                    Functional reads (in-frame and no stop codons)
``avg_clone_cdr3_num_nts``        Average clonal CDR3 length in nucleotides
``avg_clone_v_identity``          Average clonal V-region identity
``clones``                        Total number of clones
================================= ===========

Exporting Clones
----------------
In it's most basic form, the command to export clones is:

.. code-block:: bash

    $ immunedb_export PATH_TO_CONFIG clones

This will generate one file per sample each with one line per clone having the
fields below.  Note that ``intances``, ``copies``, ``avg_v_identity``, and
``top_copy_seq`` are for the clone in the context of that sample.  That is,
those fields may vary for the same clone in different samples.

================================= ===========
Field                             Description
================================= ===========
``clone_id``                      Database-wide unique clone identifier.  This
                                  number can be used to track clones across samples.
``subject``                       Subject in which the clone was found
``v_gene``                        V-gene of the clone
``j_gene``                        J-gene of the clone
``functional``                    If the clone is in-frame and contains no stop
                                  in the consensus (``T`` or ``F``)
``insertions``                    Insertions in the clone **(deprecated)**
``deletions``                     Deletions in the clone **(deprecated)**
``cdr3_nt``                       CDR3 nucleotide sequence
``cdr3_num_nts``                  CDR3 nucleotide sequence length
``cdr3_aa``                       CDR3 amino-acid sequence
``uniques``                       Unique sequences in the clone **overall**
``instances``                     Sequences instances in the clone in the
                                  associated sample
``copies``                        Copies in the clone in the associated sample
``germline``                      Clonal germline sequence
``parent_id``                     Parent ID **(deprecated)**
``avg_v_identity``                Average V-gene identity to germline
``top_copy_seq``                  Nucleotide sequence of top-copy sequence
================================= ===========

The ``--pool-on`` parameter can be used to change how data is aggregated.  By
default it takes the value ``sample`` (as described above) but it also accepts,
``subject``, or any custom metadata field(s).

For the purposes of illustration, assume we have samples with the associated
metadata below.

========    =======     =======     ======
sample      subject     tissue      subset
========    =======     =======     ======
sample1     S1          blood       naive
sample2     S1          spleen      naive
sample3     S1          spleen      mature
sample4     S3          blood       native
========    =======     =======     ======

Passing ``--pool-on subject`` will generate one file per subject with the clone
information aggregated across all samples in that subject.  Alternatively,
passing ``--pool-on tissue`` will generate one file per subject/tissue
combination.  You can pass multiple metadata fields to the ``--pool-on``
parameter as well.  For example ``--pool-on tissue subset`` will generate one
file per subject/tissue/subset combination.

Two other common parameters are ``--sample-ids`` which restricts which samples
to include in the export and ``--format`` which accepts ``immunedb`` (the
default) or ``vdjtools`` for interoperability with the `VDJtools suite
<https://vdjtools-doc.readthedocs.io>`_.

Exporting Sequences
-------------------
Sequences can be exported in `Change-O
<https://changeo.readthedocs.io/en/stable/standard.html>`_ and `AIRR
<http://docs.airr-community.org/en/latest/datarep/rearrangements.html>`_
formats.

The basic command is:

.. code-block:: bash

    $ immunedb_export PATH_TO_CONFIG sequences

This will generate one file per sample in Change-O format.  To use AIRR format,
specify ``--format airr``.  You can filter out sequences that were not
assigned to a clone with the ``--clones-only`` flag.

Exporting Selection Pressure
----------------------------
If selection pressure was calculated with the ``immunedb_clone_pressure``
command, the results can be exported in TSV format, one row per clone/sample
combination.  Additionally, unless the ``--filter samples`` is passed, there
will be one additional row per clone with a ``All Samples`` value for the
sample which indicates the overall selection pressure on the clone.

For more information on interpreting the values see `Uduman, et al, 2011
<https://www.ncbi.nlm.nih.gov/pubmed/21665923>`_ and `Yaari, et al. 2012
<https://www.ncbi.nlm.nih.gov/pubmed/22641856>`_.

========================      ==========
Field                         Value
========================      ==========
``clone_id``                  Clone ID
``subject``                   Subject to which the clone belongs
``sample``                    Sample within which the selection pressure was
                              calculated.  If ``All Samples`` the overall selection
                              pressure for the clone.
``threshold``                 The threshold at which the selection pressure was
                              calculated
``expected_REGION_TYPE``      The expected number of ``TYPE`` (``r`` or ``s``)
                              mutations in ``REGION`` (``cdr`` or ``fwr``)
``observed_REGION_TYPE``      The observed number of ``TYPE`` (``r`` or ``s``)
                              mutations in ``REGION`` (``cdr`` or ``fwr``)
``sigma_REGION``              The selection pressure in ``REGION``
``sigma_REGION_cilower``      The lower bound of the confidence interval of
                              selection in ``REGION``
``sigma_REGION_ciupper``      The upper bound of the confidence interval of
                              selection in ``REGION``
``sigma_p_REGION``            The P-value of the selection in ``REGION``
========================      ==========

Exporting MySQL Data
--------------------
The final method of exporting data is to dump the entire MySQL database to a
file.  This is meant to be a backup method rather than for downstream-analysis.

To backup run:

.. code-block:: bash

    $ immunedb_admin backup PATH_TO_CONFIG BACKUP_PATH

To restore a backup run:

.. code-block:: bash

    $ immunedb_admin restore PATH_TO_CONFIG BACKUP_PATH
