import dnautils
import json
import multiprocessing as mp
import os
import re
import traceback

from Bio import SeqIO

from sqlalchemy import distinct
from sqlalchemy.sql import desc, func

import sldb.common.config as config
import sldb.common.modification_log as mod_log
from sldb.common.models import (HashExtension, DuplicateSequence, NoResult,
                                Sample, Sequence, Study, Subject)
from sldb.identification import AlignmentException
from sldb.identification.vdj_sequence import VDJSequence
from sldb.identification.v_genes import VGermlines
from sldb.identification.j_genes import JGermlines
import sldb.util.concurrent as concurrent
import sldb.util.funcs as funcs
import sldb.util.lookups as lookups


class SampleMetadata(object):
    def __init__(self, specific_config, global_config=None):
        self._specific = specific_config
        self._global = global_config

    def get(self, key, require=True):
        if key in self._specific:

            return self._specific[key]
        elif self._global is not None and key in self._global:
            return self._global[key]
        if require:
            raise Exception(('Could not find metadata for key {}'.format(key)))


class IdentificationWorker(concurrent.Worker):
    def __init__(self, session, v_germlines, j_germlines, trim, max_vties,
                 min_similarity, sync_lock):
        self._session = session
        self._v_germlines = v_germlines
        self._j_germlines = j_germlines
        self._trim = trim
        self._min_similarity = min_similarity
        self._max_vties = max_vties
        self._sync_lock = sync_lock

    def do_task(self, args):
        meta = args['meta']
        self._print('Starting sample {}'.format(meta.get('sample_name')))
        study, sample = self._setup_sample(meta)

        vdjs = {}
        parser = SeqIO.parse(os.path.join(args['path'], args['fn']), 'fasta' if
                             args['fn'].endswith('.fasta') else 'fastq')
        # Collapse identical sequences
        self._print('\tCollapsing identical sequences')
        for record in parser:
            seq = str(record.seq)[self._trim:]
            if seq not in vdjs:
                vdjs[seq] = VDJSequence(
                    ids=[],
                    seq=seq,
                    v_germlines=self._v_germlines,
                    j_germlines=self._j_germlines,
                    quality=record.letter_annotations.get('phred_quality')
                )
            vdjs[seq].ids.append(record.description)

        self._print('\tAligning {} unique sequences'.format(len(vdjs)))
        # Attempt to align all unique sequences
        for sequence in funcs.periodic_commit(self._session, vdjs.keys()):
            vdj = vdjs[sequence]
            del vdjs[sequence]
            try:
                # The alignment was successful.  If the aligned sequence
                # already exists, append the seq_ids.  Otherwise add it as a
                # new unique sequence.
                vdj.analyze()
                if vdj.sequence in vdjs:
                    vdjs[vdj.sequence].ids += vdj.ids
                else:
                    vdjs[vdj.sequence] = vdj
            except AlignmentException:
                self.add_as_noresult(vdj, sample)
            except:
                self._print('\tUnexpected error processing sequence '
                            '{}\n\t{}'.format(vdj.ids[0],
                                              traceback.format_exc()))

        if len(vdjs) > 0:
            avg_len = sum(
                map(lambda vdj: vdj.v_length, vdjs.values())
            ) / float(len(vdjs))
            avg_mut = sum(
                map(lambda vdj: vdj.mutation_fraction, vdjs.values())
            ) / float(len(vdjs))

            self._print('\tRe-aligning {} sequences to V-ties, Mutations={}, '
                        'Length={}'.format(
                            len(vdjs), round(avg_mut, 2), round(avg_len, 2)))

            bucketed_seqs = {}
            for vdj in funcs.periodic_commit(self._session,
                                                vdjs.values()):
                del vdjs[vdj.sequence]
                try:
                    self._realign_sequence(vdj, avg_len, avg_mut)
                    bucket_key = (
                        funcs.format_ties(vdj.v_gene, 'IGHV'),
                        funcs.format_ties(vdj.j_gene, 'IGHJ'),
                        len(vdj.cdr3)
                    )
                    if bucket_key not in bucketed_seqs:
                        bucketed_seqs[bucket_key] = {}
                    bucket = bucketed_seqs[bucket_key]

                    if vdj.sequence in bucket:
                        bucket[vdj.sequence].ids += vdj.ids
                    else:
                        bucket[vdj.sequence] = vdj
                except AlignmentException:
                    self.add_as_noresult(vdj, sample)
                except:
                    self._print('\tUnexpected error processing sequence '
                                '{}\n\t{}'.format(vdj.ids[0],
                                                  traceback.format_exc()))

            # Collapse sequences that are the same except for Ns
            self._print('\tCollapsing ambiguous character sequences')
            for sequences in funcs.periodic_commit(self._session,
                                                   bucketed_seqs.values()):
                sequences = sorted(sequences.values(), cmp=lambda a, b:
                                   cmp(len(a.ids), len(b.ids)))
                while len(sequences) > 0:
                    larger = sequences.pop(0)
                    for i in reversed(range(len(sequences))):
                        smaller = sequences[i]

                        if dnautils.equal(larger.sequence, smaller.sequence):
                            larger.ids += smaller.ids
                            del sequences[i]
                    self.add_as_sequence(larger, sample, meta.get('paired'))

        sample.status = 'identified'
        self._session.commit()
        self._print('Completed sample {}'.format(sample.name))

    def cleanup(self):
        self._print('Identification worker terminating')
        self._session.close()

    def _setup_sample(self, meta):
        self._sync_lock.acquire()
        study, new = funcs.get_or_create(
            self._session, Study, name=meta.get('study_name'))

        if new:
            self._print('\tCreated new study "{}"'.format(study.name))
            self._session.commit()

        name = meta.get('sample_name')
        sample, new = funcs.get_or_create(
            self._session, Sample, name=name, study=study)
        if new:
            sample.date = meta.get('date')
            self._print('\tCreated new sample "{}"'.format(
                sample.name))
            for key in ('subset', 'tissue', 'disease', 'lab', 'experimenter',
                        'ig_class', 'v_primer', 'j_primer'):
                setattr(sample, key, meta.get(key, require=False))
            subject, new = funcs.get_or_create(
                self._session, Subject, study=study,
                identifier=meta.get('subject'))
            sample.subject = subject
            self._session.commit()

        self._sync_lock.release()

        return study, sample

    def _realign_sequence(self, vdj, avg_len, avg_mut):
        vdj.align_to_germline(avg_len, avg_mut)
        if (vdj.v_match / float(vdj.v_length) < self._min_similarity or
                len(vdj.v_gene) > self._max_vties):
            raise AlignmentException('V-match too low or too many V-ties')

    def add_as_noresult(self, vdj, sample):
        try:
            quality = funcs.ord_to_quality(vdj.quality)
            self._session.bulk_save_objects([
                NoResult(
                    seq_id=seq_id,
                    sample=sample,
                    sequence=vdj.sequence,
                    quality=quality
                ) for seq_id in vdj.ids
            ])
        except ValueError as ex:
            pass

    def add_as_sequence(self, vdj, sample, paired):
        try:
            seq = Sequence(
                seq_id=vdj.ids[0],
                sample_id=sample.id,

                subject_id=sample.subject.id,

                paired=paired,
                partial=vdj.partial,

                probable_indel_or_misalign=vdj.has_possible_indel,
                deletions=vdj.deletions,
                insertions=vdj.insertions,

                v_gene=funcs.format_ties(vdj.v_gene, 'IGHV'),
                j_gene=funcs.format_ties(vdj.j_gene, 'IGHJ'),

                num_gaps=vdj.num_gaps,
                pad_length=vdj.pad_length,

                v_match=vdj.v_match,
                v_length=vdj.v_length,
                j_match=vdj.j_match,
                j_length=vdj.j_length,

                removed_prefix=vdj.removed_prefix,
                removed_prefix_qual=funcs.ord_to_quality(
                    vdj.removed_prefix_qual),
                v_mutation_fraction=vdj.mutation_fraction,

                pre_cdr3_length=vdj.pre_cdr3_length,
                pre_cdr3_match=vdj.pre_cdr3_match,
                post_cdr3_length=vdj.post_cdr3_length,
                post_cdr3_match=vdj.post_cdr3_match,

                in_frame=vdj.in_frame,
                functional=vdj.functional,
                stop=vdj.stop,
                copy_number=len(vdj.ids),

                cdr3_nt=vdj.cdr3,
                cdr3_num_nts=len(vdj.cdr3),
                cdr3_aa=lookups.aas_from_nts(vdj.cdr3),

                sequence=str(vdj.sequence),
                quality=funcs.ord_to_quality(vdj.quality),

                germline=vdj.germline)
            self._session.add(seq)

            # Add duplicate sequences
            try:
                self._session.bulk_save_objects([
                    DuplicateSequence(
                        seq_id=seq_id,
                        duplicate_seq=seq
                    ) for seq_id in vdj.ids[1:]
                ])
            except ValueError as ex:
                pass
        except ValueError as ex:
            self.add_as_noresult(vdj, sample)


def run_identify(session, args):
    mod_log.make_mod('identification', session=session, commit=True,
                     info=vars(args))
    session.close()
    # Load the germlines from files
    v_germlines = VGermlines(args.v_germlines)
    j_germlines = JGermlines(args.j_germlines, args.upstream_of_cdr3,
                             args.anchor_len, args.min_anchor_len)
    tasks = concurrent.TaskQueue()

    sample_names = set([])
    fail = False
    for directory in args.sample_dirs:
        # If metadata is not specified, assume it is "metadata.json" in the
        # directory
        if args.metadata is None:
            meta_fn = os.path.join(directory, 'metadata.json')
        else:
            meta_fn = args.metadata

        # Verify the metadata file exists
        if not os.path.isfile(meta_fn):
            print 'Metadata file not found.'
            return

        with open(meta_fn) as fh:
            metadata = json.load(fh)

        # Create the tasks for each file
        for fn in sorted(metadata.keys()):
            if fn == 'all':
                continue
            meta = SampleMetadata(
                metadata[fn],
                metadata['all'] if 'all' in metadata else None)
            if session.query(Sample).filter(
                    Sample.name == meta.get('sample_name')
                    ).first() is not None:
                print 'Sample {} already exists. {}'.format(
                    meta.get('sample_name'), 'Skipping.' if
                    args.warn_existing else 'Cannot continue.'
                )
                fail = True
            elif meta.get('sample_name') in sample_names:
                print ('Sample {} exists more than once in metadata. Cannot '
                       'continue.').format(meta.get('sample_name'))
                return
            else:
                tasks.add_task({
                    'path': directory,
                    'fn': fn,
                    'meta': meta
                })
                sample_names.add(meta.get('sample_name'))

        if fail and not args.warn_existing:
            print ('Encountered errors.  Not running any identification.  To '
                   'skip samples that are already in the database use '
                   '--warn-existing.')
            return
    lock = mp.RLock()
    for i in range(0, min(args.nproc, tasks.num_tasks())):
        worker_session = config.init_db(args.db_config)
        tasks.add_worker(IdentificationWorker(
            worker_session, v_germlines, j_germlines, args.trim,
            args.max_vties, args.min_similarity / float(100), lock))

    tasks.start()
