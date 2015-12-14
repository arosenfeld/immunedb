# CHANGELOG
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
* `sldb_clones` now accepts the `--order` flag which will sort sequences by copy
  number for clonal assignment.

## v0.6.1.0
* Clone stats are now properly updated when `--force` flag is passed.
* Indels are now flagged for percentage mismatch in addition to windowed mutations.
* The `get_stats` API call now allows for stats to be grouped by any attribute defined in the `Sample` model.
* `sldb_sample_stats` now takes a `--clones-only` flag to only regenerate clone statistics for samples.
* CDR3 AA and CDR3 length are now properly exported.
* Lineage trees generated by neighbor joining can no longer have a zero-mutation node as the root.
* The `v_usage` API call now provides a list of groupings for the selected samples and sorts the returned V-genes.

## v0.6.0.3
* All duplicate sequences are now properly assigned an entry in the `Sequences` table, removing cycles from `DuplicateSequences`.
* Clones can now have a different V gene assigned or gaps added manually via the `sldb_modify_clone` command.
  * Modifications via `sldb_modify_clone` are recorded and can be fetched with the `modification_log` API call.  Other manual modification should make use of the `ModificationLog` model.
* Neighbor joining now properly calculates copy number.
* `sldb_clone_stats` can now be limited by clone ID.

## v0.6.0.2
* Rows summing total/unique sequences cross all samples can be included in clone
 exports.
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
