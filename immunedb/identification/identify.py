import multiprocessing as mp
import os
import traceback

from Bio import SeqIO

import immunedb.common.config as config
import immunedb.common.modification_log as mod_log
from immunedb.common.models import Sample, Study, Subject
from immunedb.identification import (add_as_noresult, add_uniques,
                                     AlignmentException)
from immunedb.identification.anchor import AnchorAligner
from immunedb.identification.metadata import parse_metadata, MetadataException
from immunedb.identification.vdj_sequence import VDJSequence
from immunedb.identification.genes import JGermlines, VGermlines
import immunedb.util.concurrent as concurrent
import immunedb.util.funcs as funcs
from immunedb.util.log import logger


class IdentificationProps(object):
    defaults = {
        'max_v_ties': 50,
        'min_similarity': .60,
        'max_padding': None,
        'trim_to': 0,
        'allow_cross_family': False,
        'max_insertions': 5,
        'max_deletions': 5,
    }

    def __init__(self, **kwargs):
        for prop, default in self.defaults.iteritems():
            setattr(self, prop, kwargs.get(prop, default))

    def valid_v_ties(self, alignment):
        return len(alignment.v_gene) <= self.max_v_ties

    def valid_min_similarity(self, alignment):
        return (alignment.v_match / float(alignment.v_length) >=
                self.min_similarity)

    def valid_padding(self, alignment):
        return (self.max_padding is None or
                alignment.seq_start <= self.max_padding)

    def valid_families(self, alignment):
        if self.allow_cross_family:
            return True

        family = None
        for gene in alignment.v_gene:
            if not family:
                family = gene.family
            elif gene.family != family:
                return False
        return True

    def valid_indels(self, alignment):
        return (len(alignment.insertions) <= self.max_insertions and
                len(alignment.deletions) <= self.max_deletions)

    def validate_cdr3(self, alignment):
        return alignment.cdr3_num_nts >= 9

    def validate(self, alignment):
        if not self.valid_min_similarity(alignment):
            raise AlignmentException(
                'V-identity too low {} < {}'.format(
                    alignment.v_match / float(alignment.v_length),
                    self.min_similarity))
        if not self.valid_v_ties(alignment):
            raise AlignmentException('Too many V-ties {} > {}'.format(
                len(alignment.v_gene), self.max_v_ties))
        if not self.valid_padding(alignment):
            raise AlignmentException('Too much padding {} (max {})'.format(
                alignment.seq_start, self.max_padding))
        if not self.valid_families(alignment):
            raise AlignmentException('Cross-family V-call')
        if not self.valid_indels(alignment):
            raise AlignmentException(
                'Too many indels insertions={} deletions={}'.format(
                    alignment.insertions, alignment.deletions))
        if not self.validate_cdr3(alignment):
            raise AlignmentException('CDR3 too short {}'.format(
                alignment.cdr3_num_nts))


class IdentificationWorker(concurrent.Worker):
    def __init__(self, session, v_germlines, j_germlines, props, sync_lock):
        self._session = session
        self._v_germlines = v_germlines
        self._j_germlines = j_germlines
        self._props = props
        self._sync_lock = sync_lock

    def do_task(self, args):
        meta = args['meta']
        self.info('Starting sample {}'.format(meta['sample_name']))
        study, sample = self._setup_sample(meta)

        vdjs = {}
        parser = SeqIO.parse(args['path'], 'fasta' if
                             args['path'].endswith('.fasta') else 'fastq')

        # Collapse identical sequences
        self.info('\tCollapsing identical sequences')
        for record in parser:
            try:
                seq = str(record.seq)
                if seq not in vdjs:
                    vdjs[seq] = VDJSequence(
                        ids=[],
                        sequence=seq,
                        quality=funcs.ord_to_quality(
                            record.letter_annotations.get('phred_quality')
                        )
                    )
                vdjs[seq].ids.append(record.description)
            except ValueError:
                continue

        alignments = {}
        aligner = AnchorAligner(self._v_germlines, self._j_germlines)
        self.info('\tAligning {} unique sequences'.format(len(vdjs)))
        # Attempt to align all unique sequences
        for sequence in funcs.periodic_commit(self._session,
                                              sorted(vdjs.keys())):
            vdj = vdjs[sequence]
            del vdjs[sequence]
            try:
                # The alignment was successful.  If the aligned sequence
                # already exists, append the seq_ids.  Otherwise add it as a
                # new unique sequence.
                alignment = aligner.get_alignment(vdj)
                seq_key = alignment.sequence.sequence
                if seq_key in alignments:
                    alignments[seq_key].sequence.ids.extend(
                        alignment.sequence.ids)
                else:
                    alignments[seq_key] = alignment
            except AlignmentException as e:
                add_as_noresult(self._session, vdj, sample, str(e))
            except Exception:
                self.error(
                    '\tUnexpected error processing sequence {}\n\t{}'.format(
                        vdj.ids[0], traceback.format_exc()))
        if len(alignments) > 0:
            avg_len = (
                sum([v.v_length for v in alignments.values()]) /
                float(len(alignments)))
            avg_mut = (
                sum([v.v_mutation_fraction for v in alignments.values()]) /
                float(len(alignments))
            )
            sample.v_ties_mutations = avg_mut
            sample.v_ties_len = avg_len

            self.info('\tRe-aligning {} sequences to V-ties, Mutations={}, '
                      'Length={}'.format(len(alignments),
                                         round(avg_mut, 2),
                                         round(avg_len, 2)))
            add_uniques(self._session, sample, alignments.values(),
                        self._props, aligner, avg_len, avg_mut)

        self._session.commit()
        self.info('Completed sample {}'.format(sample.name))

    def cleanup(self):
        self.info('Identification worker terminating')
        self._session.close()

    def _setup_sample(self, meta):
        self._sync_lock.acquire()
        self._session.commit()
        study, new = funcs.get_or_create(
            self._session, Study, name=meta['study_name'])

        if new:
            self.info('\tCreated new study "{}"'.format(study.name))
            self._session.commit()

        name = meta['sample_name']
        sample, new = funcs.get_or_create(
            self._session, Sample, name=name, study=study)
        if new:
            sample.date = meta['date']
            self.info('\tCreated new sample "{}"'.format(
                sample.name))
            for key in ('subset', 'tissue', 'disease', 'lab', 'experimenter',
                        'ig_class', 'v_primer', 'j_primer'):
                setattr(sample, key, meta.get(key, None))
            subject, new = funcs.get_or_create(
                self._session, Subject, study=study,
                identifier=meta['subject'])
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

    # If metadata is not specified, assume it is "metadata." in the
    # directory
    meta_fn = args.metadata if args.metadata else os.path.join(
        args.sample_dir, 'metadata.tsv')

    # Verify the metadata file exists
    if not os.path.isfile(meta_fn):
        logger.error('Metadata file not found.')
        return

    with open(meta_fn, 'rU') as fh:
        try:
            metadata = parse_metadata(session, fh, args.warn_existing,
                                      args.sample_dir)
        except MetadataException as ex:
            logger.error(ex.message)
            return

    # Create the tasks for each file
    for sample_name in sorted(metadata.keys()):
        tasks.add_task({
            'path': os.path.join(
                args.sample_dir, metadata[sample_name]['file_name']),
            'meta': metadata[sample_name]
        })

    props = IdentificationProps(**args.__dict__)
    lock = mp.Lock()
    for i in range(0, min(args.nproc, tasks.num_tasks())):
        worker_session = config.init_db(args.db_config)
        tasks.add_worker(IdentificationWorker(worker_session, v_germlines,
                                              j_germlines, props, lock))

    tasks.start()
