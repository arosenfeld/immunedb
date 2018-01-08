from collections import OrderedDict
import multiprocessing as mp
import os
import sys
import traceback
from sqlalchemy import func
import Queue

from Bio import SeqIO

import dnautils

import immunedb.common.config as config
import immunedb.common.modification_log as mod_log
from immunedb.common.models import (Sample, SampleMetadata, NoResult, Study,
                                    Subject)
from immunedb.identification import (add_as_noresult, add_uniques,
                                     AlignmentException)
from immunedb.identification.anchor import AnchorAligner
from immunedb.identification.metadata import (MetadataException,
                                              parse_metadata, REQUIRED_FIELDS)
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
        'genotyping': False
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


def setup_sample(session, meta):
    study, new = funcs.get_or_create(session, Study, name=meta['study_name'])

    if new:
        logger.info('Created new study "{}"'.format(study.name))
        session.commit()

    name = meta['sample_name']
    sample, new = funcs.get_or_create(session, Sample, name=name, study=study)

    if new:
        subject, new = funcs.get_or_create(
            self._session, Subject, study=study,
            identifier=meta['subject'])
        sample.subject = subject

        self.info('\tCreated new sample "{}"'.format(sample.name))
        for key, value in meta.iteritems():
            if key not in REQUIRED_FIELDS:
                self._session.add(SampleMetadata(
                    sample=sample,
                    key=key,
                    value=value
                ))


    session.commit()
    return sample


def process_vdj(vdj, aligner):
    try:
        alignment = aligner.get_alignment(vdj)
        return {
            'status': 'success',
            'alignment': alignment
        }
    except AlignmentException as e:
        return {
            'status': 'noresult',
            'vdj': vdj,
            'reason': str(e)
        }
    except Exception as e:
        return {
            'status': 'error',
            'vdj': vdj,
            'reason': str(e)
        }


def aggregate_vdj(aggregate_queue, alignments_out):
    alignments = {
        'success': {},
        'noresult': []
    }
    while True:
        try:
            result = aggregate_queue.get()
        except Queue.Empty:
            break
        if result['status'] == 'success':
            alignment = result['alignment']
            seq_key = alignment.sequence.sequence
            if seq_key in alignments['success']:
                alignments['success'][seq_key].sequence.ids.extend(
                    alignment.sequence.ids)
            else:
                alignments['success'][seq_key] = alignment
        elif result['status'] == 'noresult':
            alignments['noresult'].append(result)
        elif result['status'] == 'error':
            logger.error(
                'Unexpected error processing sequence {}\n\t{}'.format(
                    result['vdj'].ids[0], result['reason']))
    alignments['success'] = alignments['success'].values()
    alignments_out.update(alignments)


def process_vties(alignment, aligner, avg_len, avg_mut, props):
    try:
        aligner.align_to_germline(alignment, avg_len, avg_mut)
        if props.trim_to:
            alignment.trim_to(props.trim_to)

        props.validate(alignment)
        return {
            'status': 'success',
            'alignment': alignment
        }
    except AlignmentException as e:
        return {
            'status': 'noresult',
            'alignment': alignment,
            'reason': str(e)
        }
    except Exception as e:
        return {
            'status': 'error',
            'alignment': alignment,
            'reason': str(e)
        }


def aggregate_vties(aggregate_queue, seqs_out):
    bucketed_seqs = {
        'success': {},
        'noresult': []
    }
    while True:
        try:
            result = aggregate_queue.get()
        except Queue.Empty:
            break
        if result['status'] == 'success':
            alignment = result['alignment']
            bucket_key = (
                funcs.format_ties(alignment.v_gene),
                funcs.format_ties(alignment.j_gene),
                len(alignment.cdr3)
            )

            bucket = bucketed_seqs['success'].setdefault(bucket_key, {})
            if alignment.sequence.sequence in bucket:
                bucket[alignment.sequence.sequence].sequence.ids += (
                    alignment.sequence.ids
                )
            else:
                bucket[alignment.sequence.sequence] = alignment
        elif result['status'] == 'noresult':
            bucketed_seqs['noresult'].append(result)
        elif result['status'] == 'error':
            logger.error(
                'Unexpected error processing sequence {}\n\t{}'.format(
                    result['alignment'].sequence.ids[0]))

    bucketed_seqs['success'] = [
        b.values() for b in bucketed_seqs['success'].values()
    ]
    seqs_out.update(bucketed_seqs)


def process_collapse(sequences):
    sequences = sorted(
        sequences,
        key=lambda s: (len(s.sequence.ids), s.sequence.ids[0]),
        reverse=True
    )
    uniques = []
    while len(sequences) > 0:
        larger = sequences.pop(0)
        for i in reversed(range(len(sequences))):
            smaller = sequences[i]

            if dnautils.equal(larger.sequence.sequence,
                              smaller.sequence.sequence):
                larger.sequence.ids += smaller.sequence.ids
                del sequences[i]
        uniques.append(larger)
    return uniques


def aggregate_collapse(aggregate_queue, _, session, sample, props):
    while True:
        try:
            alignment = aggregate_queue.get()
            for a in alignment:
                add_as_sequence(session, a, sample,
                                strip_alleles=not props.genotyping)
        except Queue.Empty:
            break
    session.commit()


def process_sample(session, v_germlines, j_germlines, path, meta, props,
                   nproc):
    logger.info('Starting sample {}'.format(meta['sample_name']))
    sample = setup_sample(session, meta)

    vdjs = {}
    parser = SeqIO.parse(path, 'fasta' if path.endswith('.fasta') else 'fastq')

    # Collapse identical sequences
    logger.info('Collapsing identical sequences')
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

    logger.info('Aligning {} unique sequences'.format(len(vdjs)))
    aligner = AnchorAligner(v_germlines, j_germlines)

    # Initial VJ assignment
    alignments = concurrent.process_data(
        vdjs.values(),
        process_vdj,
        aggregate_vdj,
        nproc,
        process_args={'aligner': aligner},
        log_progress=True
    )
    for result in alignments['noresult']:
        add_as_noresult(session, result['vdj'], sample, result['reason'])

    alignments = alignments['success']
    if alignments:
        avg_len = (
            sum([v.v_length for v in alignments]) /
            float(len(alignments)))
        avg_mut = (
            sum([v.v_mutation_fraction for v in alignments]) /
            float(len(alignments))
        )
        sample.v_ties_mutations = avg_mut
        sample.v_ties_len = avg_len
        logger.info('Re-aligning {} sequences to V-ties: Mutations={}, '
                    'Length={}'.format(len(alignments),
                                       round(avg_mut, 2),
                                       round(avg_len, 2)))
        # Realign to V-ties
        v_ties = concurrent.process_data(
            alignments,
            process_vties,
            aggregate_vties,
            nproc,
            process_args={'aligner': aligner, 'avg_len': avg_len, 'avg_mut':
                          avg_mut, 'props': props},
            log_progress=True
        )

        for result in v_ties['noresult']:
            add_as_noresult(session, result['alignment'].sequence,
                            sample, result['reason'])

        concurrent.process_data(
            v_ties['success'],
            process_collapse,
            aggregate_collapse,
            nproc,
            aggregate_args={'session': session, 'sample': sample,
                            'props': props}
        )

        identified = int(session.query(
            func.sum(Sequence.copy_number)
        ).filter(
            Sequence.sample == sample
        ).scalar() or 0)
        noresults = int(session.query(
            func.count(NoResult.pk)
        ).filter(
            NoResult.sample == sample
        ).scalar() or 0)
        logger.info('Completed sample {} - {}/{} ({}%) identified'.format(
            sample.name,
            identified,
            identified + noresults,
            int(100 * identified / float(identified + noresults))
        ))


def run_identify(session, args):
    mod_log.make_mod('identification', session=session, commit=True,
                     info=vars(args))
    session.close()
    # Load the germlines from files
    v_germlines = VGermlines(args.v_germlines, no_ties=args.genotyping)
    j_germlines = JGermlines(args.j_germlines, args.upstream_of_cdr3,
                             args.anchor_len, args.min_anchor_len,
                             no_ties=args.genotyping)
    tasks = concurrent.TaskQueue()

    # If metadata is not specified, assume it is "metadata." in the
    # directory
    meta_fn = args.metadata if args.metadata else os.path.join(
        args.sample_dir, 'metadata.tsv')

    # Verify the metadata file exists
    if not os.path.isfile(meta_fn):
        logger.error('Metadata file not found.')
        sys.exit(-1)

    with open(meta_fn, 'rU') as fh:
        try:
            metadata = parse_metadata(session, fh, args.warn_existing,
                                      args.warn_missing, args.sample_dir)
        except MetadataException as ex:
            logger.error(ex.message)
            sys.exit(-1)

    # Create the tasks for each file
    props = IdentificationProps(**args.__dict__)
    for sample_name in sorted(metadata.keys()):
        process_sample(
            session, v_germlines, j_germlines,
            os.path.join(
                args.sample_dir,
                metadata[sample_name]['file_name']
            ),
            metadata[sample_name],
            props,
            args.nproc
        )
