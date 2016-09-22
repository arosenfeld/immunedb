import json
import multiprocessing as mp
import os
import traceback

from Bio import SeqIO

from sqlalchemy.sql import exists

import immunedb.common.config as config
import immunedb.common.modification_log as mod_log
from immunedb.common.models import Sample, Sequence, Study, Subject
from immunedb.identification import (add_as_noresult, add_uniques,
                                     AlignmentException)
from immunedb.identification.vdj_sequence import VDJSequence
from immunedb.identification.v_genes import VGermlines
from immunedb.identification.j_genes import JGermlines
import immunedb.util.concurrent as concurrent
import immunedb.util.funcs as funcs
from immunedb.util.log import logger


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
    def __init__(self, session, v_germlines, j_germlines, trim_to, max_padding,
                 max_vties, min_similarity, sync_lock):
        self._session = session
        self._v_germlines = v_germlines
        self._j_germlines = j_germlines
        self._trim_to = trim_to
        self._max_padding = max_padding
        self._min_similarity = min_similarity
        self._max_vties = max_vties
        self._sync_lock = sync_lock

    def do_task(self, args):
        meta = args['meta']
        self.info('Starting sample {}'.format(meta.get('sample_name')))
        study, sample = self._setup_sample(meta)

        vdjs = {}
        parser = SeqIO.parse(os.path.join(args['path'], args['fn']), 'fasta' if
                             args['fn'].endswith('.fasta') else 'fastq')

        # Collapse identical sequences
        self.info('\tCollapsing identical sequences')
        for record in parser:
            seq = str(record.seq)
            if seq not in vdjs:
                vdjs[seq] = VDJSequence(
                    ids=[],
                    seq=seq,
                    v_germlines=self._v_germlines,
                    j_germlines=self._j_germlines,
                    quality=funcs.ord_to_quality(
                        record.letter_annotations.get('phred_quality')
                    )
                )
            vdjs[seq].ids.append(record.description)

        self.info('\tAligning {} unique sequences'.format(len(vdjs)))
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
            except AlignmentException as e:
                add_as_noresult(self._session, vdj, sample, str(e))
            except:
                self.error(
                    '\tUnexpected error processing sequence {}\n\t{}'.format(
                        vdj.ids[0], traceback.format_exc()))
        if len(vdjs) > 0:
            avg_len = sum(
                map(lambda vdj: vdj.v_length, vdjs.values())
            ) / float(len(vdjs))
            avg_mut = sum(
                map(lambda vdj: vdj.mutation_fraction, vdjs.values())
            ) / float(len(vdjs))
            sample.v_ties_mutations = avg_mut
            sample.v_ties_len = avg_len

            self.info('\tRe-aligning {} sequences to V-ties, Mutations={}, '
                      'Length={}'.format(
                            len(vdjs), round(avg_mut, 2), round(avg_len, 2)))
            add_uniques(self._session, sample, vdjs.values(), avg_len, avg_mut,
                        self._min_similarity, self._max_vties, self._trim_to,
                        self._max_padding)

        self._session.commit()
        self.info('Completed sample {}'.format(sample.name))

    def cleanup(self):
        self.info('Identification worker terminating')
        self._session.close()

    def _setup_sample(self, meta):
        self._sync_lock.acquire()
        self._session.commit()
        study, new = funcs.get_or_create(
            self._session, Study, name=meta.get('study_name'))

        if new:
            self.info('\tCreated new study "{}"'.format(study.name))
            self._session.commit()

        name = meta.get('sample_name')
        sample, new = funcs.get_or_create(
            self._session, Sample, name=name, study=study)
        if new:
            sample.date = meta.get('date')
            self.info('\tCreated new sample "{}"'.format(
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
            logger.error('Metadata file not found.')
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
                    Sample.name == meta.get('sample_name'),
                    exists().where(
                        Sequence.sample_id == Sample.id
                    )).first() is not None:
                log_f = logger.warning if args.warn_existing else logger.error
                log_f('Sample {} already exists. {}'.format(
                    meta.get('sample_name'), 'Skipping.' if
                    args.warn_existing else 'Cannot continue.'
                ))
                fail = True
            elif meta.get('sample_name') in sample_names:
                logger.error(
                    'Sample {} exists more than once in metadata.'.format(
                        meta.get('sample_name')))
                return
            else:
                tasks.add_task({
                    'path': directory,
                    'fn': fn,
                    'meta': meta
                })
                sample_names.add(meta.get('sample_name'))

        if fail and not args.warn_existing:
            logger.error('Encountered errors.  Not running any identification.'
                         ' To skip samples that are already in the database '
                         'use --warn-existing.')
            return

    lock = mp.Lock()
    for i in range(0, min(args.nproc, tasks.num_tasks())):
        worker_session = config.init_db(args.db_config)
        tasks.add_worker(IdentificationWorker(
            worker_session, v_germlines, j_germlines, args.trim_to,
            args.max_padding, args.max_vties, args.min_similarity / float(100),
            lock))

    tasks.start()
