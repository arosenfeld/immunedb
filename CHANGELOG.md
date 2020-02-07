# CHANGELOG
## v0.29.9
* Fixed a bug where `subject` would be added to metadata after `immunedb_modify
  .. update-metadata`
* Fixed a type conversion bug for the `--reduce-difference` flag in
  `immunedb_clones`.

## v0.29.8
* Removed unused code from `immunedb_clones`.
* An IgBLAST helper script is now included in the Docker image.
* Documentation has been substantially updated and includes information about
  running IgBLAST.

## v0.29.7
* Merging similar clones in `immunedb_clones` is now much faster.

## v0.29.6
* `immunedb_identify` now accepts gzipped files.

## v0.29.5
* Minor fix for sequences with indels in `immunedb_modify ... combine-samples`

## v0.29.3
* Collapsing sequences imported from AIRR format now uses substantially less
  memory.

## v0.29.2
* CDR3 indels are now supported for AIRR formatted input.
* Clonal assignment has been greatly enhanced by allowing sequences with
  different V-gene and J-gene calls but similar CDR3 nucleotide sequences to be
  combined.  This is helpful in highly mutated populations where gene
  assignments may be incorrect for some sequences.
* Clones with less than or equal to 4 (default) nucleotide differences are
  automatically combined.  You can tweak this behavior in `immunedb_clones`
  with the `--reduce-difference` flag.

## v0.29.0
* `immunedb_import` for IgBLAST results in AIRR format has been greatly
  improved both in performance and accuracy.
* Short-read edge case for J-alignment now assigns proper gene.
* `immunedb_identify` now requires the `--ties` flag to be specified to enable
  V-tie calculations.  By default simply the highest-identity gene is assigned.
* `immunedb_clone_trees` defaults to only include the V-gene in lineage
  construction.  To use full sequences, specify the `--full-seq` flag.
* `immunedb_clones` now implements hierarchical clustering with the `cluster`
  method.
* The `lineage` method for `immunedb_clones` has been removed.
* `immunedb_clones` now defaults to resetting existing clones.  To avoid this
  behavior specify the `--skip-regen` flag.
* `immunedb_clones` now automatically combines clones with the same V, J, and
  CDR3 AA sequence.  To avoid this behavior specify the `--skip-reduce` flag.
* `immunedb_clones` now defaults to assigning sub-clones.

## v0.28.3
* Importing is now supported for only IgBLAST in AIRR format.
* Database configuration can now be specified with environment variables.
  `IMMUNEDB_DB`, `IMMUNEDB_PASS`, `IMMUNEDB_USER` are required to be set and
  `IMMUNEDB_HOST` is optional (defaults to `localhost`).  When using
  environment variables, the common first parameter to all commands,
  `db_config` should be excluded.
* Sample export now includes average clonal V-identity and CDR3 length as
  columns.

## v0.28.2
* J-assignment fix for reverse-complement sequences.

## v0.28.1
* Minor fix for V-identity when exporting clones.

## v0.28.0
* Fix for clonal similarity.
* The Docker container now proxies both the frontend and API through port 8080.
  Instead of running `serve_immunedb.sh`, simply navigate to
  http://localhost:8080/frontend/`db_name`.
* Exporting clones now includes a `v_identity` column.

## v0.27.0
* Identity filtering for all CLI commands has been changed to use fractions
  rather than integers.
* The `similarity` clonal assignment method now has a `--level` flag taking
  either `aa` (the default) or `nt` which sets the level at which similarity is
  calculated.
* Due to the change above, the `tcell` clonal assignment method has been
  removed.  Use `similarity` with the flags `--level nt --min-similarity 1`
  for the same results.

## v0.26.0
* J-genes are now assigned based on percent identity rather than number of
  matching bases which favored longer J-genes.
* Exporting has been cleaned up and the API and CLI have been unified.
* Clones can now be pooled when exported.
* Exporting via the API now returns a `uid` token rather than the result.  The
  result can be accessed once the data is ready using the `uid`.
* A new `overlap` export format has been added both to the CLI and
  web-interface.

## v0.25.0
* Exporting has been entirely rewritten.  Exporting from the CLI and
  web-interface via the API now use unified methods.
* Duplicate sequences at the sample level are no longer stored.  This increase
  identification / importing speed drastically.
* Selection pressure, sample metadata, clones summaries, and clonal overlap
  have all been added to `immunedb_export`.
* `immunedb_clone_trees` now has new parameters to filter mutations & sequences
  by copy number and number of samples.
* Importing for Adaptive data has now been added to the `immunedb_import`
  command.
* Metadata can now be updated and samples can be combined with the
  `immunedb_modify` command.

## v0.24.0
* ImmuneDB has now been ported to Python 3.  ImmuneDB v0.23.0 will be
  maintained in a separate branch for bug fixes only.  No new Python 2 features
  will be added.
* Multiprocessing for identification has been rewritten to process each sample
  sequentially.
* Fixed a bug where clone counts were incorrect during API call.
* Adding option to export in AIRR format.
* Clearcut must now be located in $PATH.
* The Docker container has been rewritten and no longer uses Docker Compose.

## v0.23.0
* Arbitrary metadata can now be specified for samples.  The only required
  fields are `sample_name`, `subject`, and `study_name`.
* Identification now uses multiprocessing for each sample independently, rather
  than one process per sample.
* Sequences with stop codons in the CDR3 can now be added to clones.
* V-length is now correctly calculated for 5' trimmed sequences.
* `immunedb_sql` now has the optional `--query` argument to run a query via the
  CLI.
* Exporting via the CLI should now be faster.
* Insertions and deletions can now be exported for clones and sequences.
* Selection pressure is now provided with the sample analysis API call.

## v0.22.0
* Arbitary fields can be specified in metadata files.  The only required fields
  are now `file_name`, `study_name`, `sample_name`, and `subject`.
* Clones correctly filtered when functionality is not stipulated.
* Sequences that comprise trees can now be filtered by their copy number using
  `--min-seq-copies`.
* Sequences with invalid bases will no longer cause identification to fail.
* Various local-alignment bug fixes.
* The number of instances in a clone is now stored in the database.
* Mutations are correctly exported when limiting to certain samples.
* Sequences can now be exported in AIRR compatible GenBank format with
  `immunedb_export genbank`.
* Sequences can now be exported in Change-O format with `immunedb_export
  changeo`.
* Genotyping with TIgGER can be run by passing the `--genotype` flag to
  `immunedb_identify`.  See the docs for more instructions.
* Selection pressure is now stored in a separate table for easier querying.
* Optional baseline regression tests added.
* Tree node features with no values are now excluded.

## v0.21.0
* Local alignment has been entirely rewritten to use bowtie2.  This drastically
  reduces the time necessary to locally align sequences.
* The `/samples/overlap/` API endpoint now properly returns clones when not
  filtering by functionality.
* A new `--min-seq-copies` flag has been added to `immunedb_clone_trees` which
  limits which sequences are included in trees based on their copy number.
* Re-creating trees of specifically specified clones now requires the `--force`
  flag.
* Default arguments for all commands are now automatically populated.
* API calls to list subjects now provides unique sequences as well as copies
  and instances.
* Arguments for `immunedb_clones` have been changed.  Three methods are now
  available, similarity based for B-cells, identical based for T-cells, and a
  lineage-separation method.
* Instance counts are now stored for clones.

## v0.20.1
* Metadata is now specified in TSV files rather than JSON.  A
  `immunedb_metadata` command has been added to automatically create a template
  file from FASTA/FASTQ files.
* Gene representation has been modified to be more flexible with non-standard
  naming schemes.
* J-gene identification has been substantially changed.  If anchoring does not
  work, full-sequence matching is attempted.
* Clones can now be exported in the format expected by VDJtools with the
  `immunedb_export vdjtools` command.

## v0.19.1
* Subclones are now properly assigned.
* The order of exported sequences is now consistent.
* Updates for consistency with Python 3 test cases.

## v0.19.0
* ImmuneDB can now process T-cell sequences.
* Clonal assignment now includes an optional "subclone" process for
  locally-aligned sequences.  Subclones are clones which share features of their
  parent, but contain insertions or deletions.  See the documentation for more
  information.
* External clonal assignments can now be imported with `immunedb_clone_import`.
* Optional rollbar support has been added to `immunedb_rest` to track errors.
* Logging has been overhauled and is more consistent with best-practices.

## v0.18.0
* The package has been renamed from AIRRDB to ImmuneDB.
* Tests now check for the presence of a local-alignment binary.
* SciPy has been removed as a dependency and a custom hypergeom function has
  been included.

## v0.17.0
* The package has been renamed from SLDB to AIRRDB.
* J-gene offsets are now set to human values by default.
* Local alignment has been updated and should properly work for most sequences.

## v0.16.2
* Selection pressure can now be calculated for mutations happening exactly a
  specified number of times.
* Clonal overlap calculations are now faster.
* A `sldb_sql` command has been added to ease direct interface with MySQL.
* API call for clonal overlap now properly pages.
* Improved error handling in lineage construction.
* Docker compose is now used to separate the different AIRRDB components.

## v0.16.1
* J-genes are now properly assigned.

## v0.16.0
* Alleles are no longer annotated.
* Sequences can be optionally trimmed during identification or importing.
* Sequences with stop codons can optionally be excluded from lineages.
* Sequences are now properly assigned to clones regardless of CDR3 length.
* Sequences with ambiguous bases in their J-genes are now properly identified.
* V- and J-gene tie code has been consolidated.
* Total clone copy number is now properly calculated for statistics.
* Sequences with ambiguous CDR3s are now properly added to clones.
* V-ties for locally-aligned sequences are now properly annotated.
* Mutation rate for each sample is now stored in the underlying database.
* Lineage node copy numbers are now correct for collapsed sequences
* Clone overlap queries are now much faster.

## v0.15.0
* Local alignment now uses external libraries.
* Insertions and deletions are now included in sequence records
* A Dockerfile is now available for AIRRDB.

## v0.14.1
* Exporting clones by sample now works properly.
* Memory usage and run-time for local-alignment has been reduced.
* Documentation has been cleaned up.
* Selection pressure can now be calculated at any level.

## v0.14.0
* The API has been simplified and re-organized.
* URLs for API calls now use run length encoding to specify which samples to
  analyze.  This fixes issues when many samples are selected and cause the URL
  to be too long.
* Grouped quality scores are now properly calculated.
* Rarefaction calculations have been removed.

## v0.13.0
* Clone mutations can now have arbitrary thresholds.
* The clone overlap query has been optimized and now properly filters
  functional and non-functional clones.

## v0.12.0
* `sldb_admin` has been added to simplify creating, deleting, backing up, and
  restoring SLDB instances.
* `sldb_local_align` has been added for locally aligning sequences marked as
  having insertions or deletions.
* `sldb_clone_selection_pressure` has been added to calculate clonal selection
  pressure.  `sldb_clone_stats` now only calculates mutations and overlap, but
  much more quickly.
* Duplicate sequences, regardless of ambiguous characters, are automatically
  collapsed during identification.
* API call `get_stats` now takes a `percentages` parameter which will return
  statistics as percentages.
* SLDB no longer uses two databases and now only requires one configuration
  file.
* Identification speed has been increased..
* Samples can now be annotated with an `ig_class` specifying the isotype of the
  sample (e.g. IgA, IgE).
* Sequences instances are now counted at the subject level.

## v0.11.4
* Exporting clone overlap now includes selected and all samples.
* Identification will no longer fail for samples with zero identifiable reads.
* It is no longer possible to have multiple input files for one sample.
  Additionally, identification will not allow sequences to be added to existing
  samples.
* Hypergeometric probabilities for V-ties are now cached, greatly improving
  identification performance.
* A `--trim INT`` parameter has been added to identification allowing reads to
  be trimmed prior to identification.

## v0.10.0
* TokuDB has been dropped in favor of InnoDB for the purpose of easier
  installation.
* Identification tests have been re-written.
* Sequences that cannot be inserted due to a field-length restriction are added
  as `NoResult`s whenever possible.
* Clonal assignment now includes partial reads by default.
* Collapsing of sequences now occurs within V, J, CDR3 length buckets for
  efficiency.
* Identification now looks for `D.....C` in sequences if all other anchors
  fail.
* Sample-level duplicate sequences now have the correct clone ID after.
* V-match percentage is now correct for partial sequences.

## v0.9.0
* Quality strings are now properly oriented for reverse-complement sequences.
* Trees will no longer have zero-mutation roots.
* Clones can now be created with an specifiable minimum-copy number.
* Sequence exports can now optionally only include sequences assigned to clones.
* Total sequence counts in sample statistics now work with sample-level
  collapsing.
* Tree creation will now emit a warning when mutation information is
  unavailable.
* Multiprocess workers now emit warning when uncaught errors occur.
* VDJ alignment now uses exceptions to indicate alignment failures.
* String-fields in models are now verified to be of correct length or a
  `ValueEror` is thrown.
* CDR3s are now limited to the lesser of 32 amino acids or 96 nucleotides.
* Models now consistently use `cdr3` instead of `junction` for the CDR3 region.
* Identification has been refactored to be cleaner and more efficient.
* Regression testing has been added in the `tests` directory.
* The `CloneGroup` model has been removed.
* Exporting has been refactored.

## v0.8.0
* J gene germlines are now specified by a FASTA file than hard-coded sequences.
* Clone lineages can now be created only from mutations that occur in a given
  number of samples.
* Various performance enhancements to clone statistics.
* Selection pressure is pre-calculated for all mutations as well as those which
  occur at least twice.
* Removed clone collapse level since it will never result in further collapsing
  past the subject level.
* Delimited importing update to match new models.
* Sequences with various capitalization is now normalized.
* V identification no longer looks for hard-coded anchors.

## v0.7.0
* Versioning will now follow the [Semantic Versioning
  Standard](http://semver.org).
* Major Feature: Sequences must now be collapsed at three different levels: the
  sample, subject, and clone.  This collapsing is to take into account Ns added
  from quality filtering.  Sequences that are identical except for Ns are
  considered the same and will be collapsed into the highest copy-number
  sequence.  Equality checking ignoring Ns is written in C for efficiency.
* Feature: Fully aligned sequences with V and J assignments can now be imported
  from CSV files.
* Feature: Phred quality scores can now be analyzed from FASTQ files.
* Feature: Rarefaction for samples can now be calculated with the `rarefaction`
  API call.
* Feature: V-gene diversity for samples can now be calculate with the
  `diversity` API call.
* Feature: Sequences can now be discarded based on number of V-ties and a
  minimum identity-to-germline threshold during identification.
* Enhancement: When checking if a sequence has a similar CDR3 for clonal
  assignment, only unique CDR3 amino-acid sequences are checked (#24).
* Enhancement: Clone overlap in a sample-context can now be exported as a CSV
  via the `clone_overlap` API call.
* Bug fix: Workers will no longer prematurely terminate due to blocking on the
  task queue.
* Bug fix: Grouping of sample statistics no longer inflates distribution values.
* Bug fix: Copy numbers for duplicate sequences during `sldb_identify` are now
  correct.

## v0.6.2.0
* [Baseline](http://selection.med.yale.edu/baseline) has now been integrated to
  calculate clonal selection pressure during clone statistic calculations.
* Clone comparison now only allows one clone to be selected.
* Modification log messages now added at each pipeline stage.
* Mutations are now precalculated for all sequences and clones.
* Mutations can now be exported for both clones and samples.
* All pipeline stages now use the multiprocessing module to parallelize
  processing.

## v0.6.1.2
* Mutations can now be filtered by occurrence frequency via the REST API.

## v0.6.1.1
* `sldb_sample_stats` now accepts the `--clones-only` flag which, when set, will
  cause sample statistics only to be generated for clone filters.  Useful for
  updates to clonal assignment methods.
* Fixed a bug where exporting sequences did not return the CDR3 NTs, AAs, or
  length.
* V-gene names are now lexicographically sorted when requested via the REST API.
* `sldb_clones` now accepts the `--order` flag which will sort sequences by
  copy number for clonal assignment.

## v0.6.1.0
* Clone stats are now properly updated when `--force` flag is passed.
* Indels are now flagged for percentage mismatch in addition to windowed
  mutations.
* The `get_stats` API call now allows for stats to be grouped by any attribute
  defined in the `Sample` model.
* `sldb_sample_stats` now takes a `--clones-only` flag to only regenerate clone
  statistics for samples.
* CDR3 AA and CDR3 length are now properly exported.
* Lineage trees generated by neighbor joining can no longer have a
  zero-mutation node as the root.
* The `v_usage` API call now provides a list of groupings for the selected
  samples and sorts the returned V-genes.

## v0.6.0.3
* All duplicate sequences are now properly assigned an entry in the `Sequences`
  table, removing cycles from `DuplicateSequences`.
* Clones can now have a different V gene assigned or gaps added manually via
  the `sldb_modify_clone` command.
  * Modifications via `sldb_modify_clone` are recorded and can be fetched with
    the `modification_log` API call.  Other manual modification should make use
    of the `ModificationLog` model.
* Neighbor joining now properly calculates copy number.
* `sldb_clone_stats` can now be limited by clone ID.

## v0.6.0.2
* Rows summing total/unique sequences cross all samples can be included in
  clone exports.
* Neighbor joining now added as a method of lineage tree creation.
* Fixed incorrect unique sequence count in clone comparison.

## v0.6.0.0
* V gene usage API tweak to allow for exporting via website.
* Mutation frequency is now calculated for various thresholds.
* HighV-Quest output can now be imported with the `sldb_hvquest` binary.
 * Optionally, V-ties can be calculated in addition to the Vs specified
* Sequences with probable indels or misalignments are excluded from clonal
  assignment by default.  Override this behavior in `sldb_clones` with
  `--include-indels`.
* Major data-model changes:
 * Consolidated `SequenceMapping` to `Sequence` model to reduce joining.
 * `SampleStatistics` changed to match updated `Sequence` model.
 * Added `CloneStats` model to reduce API query time.
* Added binary `sldb_clone_stats` to populate `CloneStats` models.
* Renamed `aggregation/stats.py` to `aggregation/sample_stats.py` and
`sldb_stats` binary to `sldb_sample_stats` for new clone statistics scripts.
* Can now export both sequences and clones.  Re-architected exporting classes.

## v0.5.0.1
* Fixed case where duplicate sequences were incorrectly inserted.

## v0.5.0.0
* **First stable release**
* Duplicate sequences are now detected during identification and not
re-identified.
* Metadata fallthrough to "all" block properly during identification.
* pRESTO references removed in favor of "R1+R2"
* Insertion/deletion check based on sliding window.
* V and J ties now calculated on a per-sample basis.
* Major changes and fixes to V/J identification:
  * V gene alleles can now be identified and must be separated with an asterisk
    (e.g.  IGHV4-34*01).
  * Anchors are now found using reversed frame-shifting if forward
  * frame-shifting yields a no-result.
  * V and J genes now match into the CDR3 based on sliding window.
  * Germlines can now be specified during identification.  Note each germline name
    must refer to a unique sequence.
