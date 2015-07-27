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


class SequenceRecord(object):
    def __init__(self, sequence, quality):
        self.sequence = sequence
        self.quality = quality

        self.seq_ids = []
        self.vdj = None

    def add_as_noresult(self, session, sample):
        try:
            for seq_id in self.seq_ids:
                session.add(NoResult(
                    seq_id=seq_id,
                    sample=sample,
                    sequence=self.sequence))
        except ValueError as ex:
            pass

    def add_as_duplicate(self, session, sample, existing_seq_id):
        try:
            for seq_id in self.seq_ids:
                session.add(DuplicateSequence(
                    seq_id=seq_id,
                    sample=sample,
                    duplicate_seq_id=existing_seq_id))
        except ValueError as ex:
            pass

    def add_as_sequence(self, session, sample, meta):
        existing = session.query(Sequence).filter(
            Sequence.sample_seq_hash == HashExtension.hash_fields(
                (sample.id, self.vdj.sequence))
        ).first()
        if existing is not None:
            existing.copy_number += len(self.seq_ids)
            self.add_as_duplicate(session, sample, existing.seq_id)
            return

        if self.vdj.quality is not None:
            # Converts quality array into Sanger FASTQ quality string
            quality = ''.join(map(
                lambda q: ' ' if q is None else chr(q + 33), self.vdj.quality))
        else:
            quality = None

        try:
            session.add(Sequence(
                seq_id=self.vdj.id,
                sample_id=sample.id,

                paired=meta.get('paired'),
                partial=self.vdj.partial,
                probable_indel_or_misalign=self.vdj.has_possible_indel,

                v_gene=funcs.format_ties(self.vdj.v_gene, 'IGHV'),
                j_gene=funcs.format_ties(self.vdj.j_gene, 'IGHJ'),

                num_gaps=self.vdj.num_gaps,
                pad_length=self.vdj.pad_length,

                v_match=self.vdj.v_match,
                v_length=self.vdj.v_length,
                j_match=self.vdj.j_match,
                j_length=self.vdj.j_length,

                pre_cdr3_length=self.vdj.pre_cdr3_length,
                pre_cdr3_match=self.vdj.pre_cdr3_match,
                post_cdr3_length=self.vdj.post_cdr3_length,
                post_cdr3_match=self.vdj.post_cdr3_match,

                in_frame=self.vdj.in_frame,
                functional=self.vdj.functional,
                stop=self.vdj.stop,
                copy_number=len(self.seq_ids),

                cdr3_nt=self.vdj.cdr3,
                cdr3_num_nts = len(self.vdj.cdr3),
                cdr3_aa=lookups.aas_from_nts(self.vdj.cdr3),

                sequence=str(self.vdj.sequence),
                quality=quality,

                germline=self.vdj.germline))
        except ValueError as ex:
            self.add_as_noresult(session, sample)


class IdentificationWorker(concurrent.Worker):
    def __init__(self, session, v_germlines, j_germlines, max_vties,
            min_similarity, sync_lock):
        self._session = session
        self._v_germlines = v_germlines
        self._j_germlines = j_germlines
        self._min_similarity = min_similarity
        self._max_vties = max_vties
        self._sync_lock = sync_lock

    def do_task(self, args):
        meta = args['meta']
        self._print('Starting sample {}'.format(meta.get('sample_name')))
        study, sample = self._setup_sample(meta)

        sequences = {}
        parser = SeqIO.parse(os.path.join(args['path'], args['fn']), 'fasta' if
                             args['fn'].endswith('.fasta') else 'fastq')
        # Collapse identical sequences
        self._print('\tCollapsing identical sequences')
        for record in parser:
            seq = str(record.seq)
            if seq not in sequences:
                sequences[seq] = SequenceRecord(
                    seq, record.letter_annotations.get('phred_quality'))
            sequences[seq].seq_ids.append(record.description)

        self._print('\tAligning unique sequences')
        # Attempt to align all unique sequences
        for sequence in funcs.periodic_commit(self._session, sequences.keys()):
            record = sequences[sequence]
            del sequences[sequence]

            try:
                vdj = VDJSequence(
                    record.seq_ids[0],
                    record.sequence,
                    self._v_germlines,
                    self._j_germlines,
                    quality=record.quality
                )
                # The alignment was successful.  If the aligned sequence
                # already exists, append the seq_ids.  Otherwise add it as a
                # new unique sequence.
                record.vdj = vdj
                if vdj.sequence in sequences:
                    sequences[vdj.sequence].seq_ids += record.seq_ids
                else:
                    sequences[vdj.sequence] = record
            except AlignmentException:
                record.add_as_noresult(self._session, sample)
            except:
                self._print('\tUnexpected error processing sequence '
                            '{}\n\t{}'.format(record.seq_ids[0],
                                            traceback.format_exc()))

        avg_len = sum(
            map(lambda r: r.vdj.v_length, sequences.values())
        ) / float(len(sequences))
        avg_mut = sum(
            map(lambda r: r.vdj.mutation_fraction, sequences.values())
        ) / float(len(sequences))

        self._print('\tRe-aligning to V-ties, Mutations={}, Length={}'.format(
            round(avg_mut, 2), round(avg_len, 2)))
        for record in funcs.periodic_commit(self._session, sequences.values()):
            try:
                self._realign_sequence(record.vdj, avg_len, avg_mut)
                record.add_as_sequence(self._session, sample, meta)
            except AlignmentException:
                record.add_as_noresult(self._session, sample)
            except:
                self._print('\tUnexpected error processing sequence '
                            '{}\n\t{}'.format(record.seq_ids[0],
                                            traceback.format_exc()))
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
            self._print('\tCreated new sample "{}" in MASTER'.format(
                sample.name))
            for key in ('subset', 'tissue', 'disease', 'lab', 'experimenter'):
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
        if (vdj.v_match / float(vdj.v_length) < self._min_similarity
                or len(vdj.v_gene) > self._max_vties):
            raise AlignmentException('V-match too low or too many V-ties')


def run_identify(session, args):
    mod_log.make_mod('identification', session=session, commit=True,
                     info=vars(args))
    session.close()
    # Load the germlines from files
    v_germlines = VGermlines(args.v_germlines)
    j_germlines = JGermlines(args.j_germlines, args.upstream_of_cdr3,
                             args.anchor_len, args.min_anchor_len)
    tasks = concurrent.TaskQueue()

    # If metadata is not specified, assume it is "metadata.json" in the
    # samples_dir directory
    if args.metadata is None:
        meta_fn = '{}/metadata.json'.format(args.samples_dir)
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
        tasks.add_task({
            'path': args.samples_dir,
            'fn': fn,
            'meta': SampleMetadata(
                metadata[fn],
                metadata['all'] if 'all' in metadata else None
            )
        })

    lock = mp.RLock()
    for i in range(0, args.nproc):
        worker_session = config.init_db(args.master_db_config,
                                        args.data_db_config)
        tasks.add_worker(IdentificationWorker(worker_session,
                                              v_germlines,
                                              j_germlines,
                                              args.max_vties,
                                              args.min_similarity / float(100),
                                              lock))

    tasks.start()
