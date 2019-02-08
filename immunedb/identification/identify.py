import os
import sys
import time
from sqlalchemy import func

from Bio import SeqIO

import dnautils

import immunedb.common.config as config
import immunedb.common.modification_log as mod_log
from immunedb.common.models import (Sample, SampleMetadata, Sequence, NoResult,
                                    Study, Subject)
from immunedb.identification import (add_noresults_for_vdj, add_sequences,
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
        'genotyping': False,
    }

    def __init__(self, **kwargs):
        for prop, default in self.defaults.items():
            setattr(self, prop, kwargs.get(prop, default))

    def valid_v_ties(self, alignment):
        return len(alignment.v_gene) <= self.max_v_ties

    def valid_min_similarity(self, alignment):
        return alignment.v_match / alignment.v_length >= self.min_similarity

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
                    alignment.v_match / alignment.v_length,
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
            session, Subject, study=study,
            identifier=meta['subject'])
        sample.subject = subject

        for key, value in meta.items():
            if key not in REQUIRED_FIELDS:
                session.add(SampleMetadata(
                    sample=sample, key=key, value=value
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


def aggregate_vdj(aggregate_queue):
    alignments = {
        'success': {},
        'noresult': []
    }
    for result in aggregate_queue:
        if result['status'] == 'success':
            alignment = result['alignment']
            seq_key = alignment.sequence.sequence
            if seq_key in alignments['success']:
                alignments['success'][seq_key].sequence.copy_number += (
                    alignment.sequence.copy_number)
            else:
                alignments['success'][seq_key] = alignment
        elif result['status'] == 'noresult':
            alignments['noresult'].append(result)
        elif result['status'] == 'error':
            logger.error(
                'Unexpected error processing sequence {}\n\t{}'.format(
                    result['vdj'].seq_id, result['reason']))
    alignments['success'] = alignments['success'].values()
    return alignments


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


def aggregate_vties(aggregate_queue):
    bucketed_seqs = {
        'success': {},
        'noresult': []
    }
    for result in aggregate_queue:
        if result['status'] == 'success':
            alignment = result['alignment']
            bucket_key = (
                funcs.format_ties(alignment.v_gene),
                funcs.format_ties(alignment.j_gene),
                len(alignment.cdr3)
            )

            bucket = bucketed_seqs['success'].setdefault(bucket_key, {})
            if alignment.sequence.sequence in bucket:
                bucket[alignment.sequence.sequence].sequence.copy_number += (
                    alignment.sequence.copy_number
                )
            else:
                bucket[alignment.sequence.sequence] = alignment
        elif result['status'] == 'noresult':
            bucketed_seqs['noresult'].append(result)
        elif result['status'] == 'error':
            logger.error(
                'Unexpected error processing sequence {}\n\t{}'.format(
                    result['alignment'].sequence.seq_id))

    bucketed_seqs['success'] = [
        b.values() for b in bucketed_seqs['success'].values()
    ]
    return bucketed_seqs


def process_collapse(sequences):
    sequences = sorted(
        sequences,
        key=lambda s: (s.sequence.copy_number, s.sequence.seq_id),
        reverse=True
    )
    uniques = []
    while len(sequences) > 0:
        larger = sequences.pop(0)
        for i in reversed(range(len(sequences))):
            smaller = sequences[i]
            if dnautils.equal(larger.sequence.sequence,
                              smaller.sequence.sequence):
                larger.sequence.copy_number += smaller.sequence.copy_number
                del sequences[i]
        uniques.append(larger)
    return uniques


def aggregate_collapse(aggregate_queue, db_config, sample_id, props):
    seqs_to_add = []
    session = config.init_db(db_config, create=False)
    sample = session.query(Sample).filter(Sample.id == sample_id).one()
    for i, alignment in enumerate(aggregate_queue):
        for seq in alignment:
            seqs_to_add.append(seq)
            if len(seqs_to_add) >= 1000:
                add_sequences(session, seqs_to_add, sample,
                              strip_alleles=not props.genotyping)
                seqs_to_add = []
                session.commit()
    if seqs_to_add:
        add_sequences(session, seqs_to_add, sample,
                      strip_alleles=not props.genotyping)
    logger.info('Finished aggregating sequences')
    session.commit()
    session.close()


def read_input(path):
    vdjs = []
    parser = SeqIO.parse(path, 'fasta' if path.endswith('.fasta') else 'fastq')

    # Collapse identical sequences
    logger.info('Parsing input')
    for record in parser:
        try:
            vdjs.append(VDJSequence(
                seq_id=record.description,
                sequence=str(record.seq),
                quality=funcs.ord_to_quality(
                    record.letter_annotations.get('phred_quality')
                )
            ))
        except ValueError:
            continue

    logger.info('There are {} sequences'.format(len(vdjs)))
    return vdjs


def process_sample(db_config, v_germlines, j_germlines, path, meta, props,
                   nproc):
    session = config.init_db(db_config)
    start = time.time()
    logger.info('Starting sample {}'.format(meta['sample_name']))
    sample = setup_sample(session, meta)

    aligner = AnchorAligner(v_germlines, j_germlines)

    # Initial VJ assignment
    alignments = concurrent.process_data(
        read_input,
        process_vdj,
        aggregate_vdj,
        nproc,
        process_args={'aligner': aligner},
        generate_args={'path': path},
    )
    logger.info('Adding noresults')
    for result in alignments['noresult']:
        add_noresults_for_vdj(session, result['vdj'], sample, result['reason'])

    alignments = alignments['success']
    if alignments:
        avg_len = (
            sum([v.v_length for v in alignments]) /
            len(alignments))
        avg_mut = (
            sum([v.v_mutation_fraction for v in alignments]) /
            len(alignments)
        )
        sample.v_ties_mutations = avg_mut
        sample.v_ties_len = avg_len
        logger.info('Re-aligning {} sequences to V-ties: Mutations={}, '
                    'Length={}'.format(len(alignments),
                                       round(avg_mut, 2),
                                       round(avg_len, 2)))
        session.commit()
        # Realign to V-ties
        v_ties = concurrent.process_data(
            alignments,
            process_vties,
            aggregate_vties,
            nproc,
            process_args={'aligner': aligner, 'avg_len': avg_len, 'avg_mut':
                          avg_mut, 'props': props},
        )
        logger.info('Adding noresults')

        for result in funcs.periodic_commit(session, v_ties['noresult'], 100):
            add_noresults_for_vdj(session, result['alignment'].sequence,
                                  sample, result['reason'])

        logger.info('Collapsing {} buckets'.format(len(v_ties['success'])))
        session.commit()

        # TODO: Change this so we arent copying everything between processes
        concurrent.process_data(
            [list(v) for v in v_ties['success']],
            process_collapse,
            aggregate_collapse,
            nproc,
            aggregate_args={'db_config': db_config, 'sample_id': sample.id,
                            'props': props}
        )
        session.expire_all()
        session.commit()

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
        if identified + noresults:
            frac = int(100 * identified / (identified + noresults))
        else:
            frac = 0
        logger.info(
            'Completed sample {} in {}m - {}/{} ({}%) identified'.format(
                sample.name,
                round((time.time() - start) / 60., 1),
                identified,
                identified + noresults,
                frac
            )
        )
    session.close()


def run_identify(session, args):
    mod_log.make_mod('identification', session=session, commit=True,
                     info=vars(args))
    # Load the germlines from files
    v_germlines = VGermlines(args.v_germlines, no_ties=args.genotyping)
    j_germlines = JGermlines(args.j_germlines, args.upstream_of_cdr3,
                             args.anchor_len, args.min_anchor_len,
                             no_ties=args.genotyping)

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
            logger.error(ex)
            sys.exit(-1)

    session.close()
    # Create the tasks for each file
    props = IdentificationProps(**args.__dict__)
    for sample_name in sorted(metadata.keys()):
        process_sample(
            args.db_config, v_germlines, j_germlines,
            os.path.join(
                args.sample_dir,
                metadata[sample_name]['file_name']
            ),
            metadata[sample_name],
            props,
            args.nproc
        )
