# CHANGELOG

## v6.0.0
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

## v5.0.1
* Fixed case where duplicate sequences were incorrectly inserted.

## v5.0.0
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
