import json
import multiprocessing as mp
import os
import Queue
import traceback

from Bio import SeqIO

import sldb.common.config as config
import sldb.common.modification_log as mod_log
from sldb.common.models import Sample, Study, Subject
from sldb.identification import AlignmentException, SequenceRecord
from sldb.identification.vdj_sequence import VDJSequence
from sldb.identification.v_genes import VGermlines
from sldb.identification.j_genes import JGermlines
import sldb.util.concurrent as concurrent
import sldb.util.funcs as funcs


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

        sequences = {}
        parser = SeqIO.parse(os.path.join(args['path'], args['fn']), 'fasta' if
                             args['fn'].endswith('.fasta') else 'fastq')
        # Collapse identical sequences
        self._print('\tCollapsing identical sequences')
        for record in parser:
            seq = str(record.seq)[self._trim:]
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
                if record.sequence in sequences:
                    sequences[record.sequence].seq_ids += record.seq_ids
                else:
                    sequences[record.sequence] = record
            except AlignmentException:
                record.add_as_noresult(self._session, sample)
            except:
                self._print('\tUnexpected error processing sequence '
                            '{}\n\t{}'.format(record.seq_ids[0],
                                              traceback.format_exc()))

        if len(sequences) > 0:
            avg_len = sum(
                map(lambda r: r.vdj.v_length, sequences.values())
            ) / float(len(sequences))
            avg_mut = sum(
                map(lambda r: r.vdj.mutation_fraction, sequences.values())
            ) / float(len(sequences))

            self._print('\tRe-aligning to V-ties, Mutations={}, '
                        'Length={}'.format(
                            round(avg_mut, 2), round(avg_len, 2)))

            for record in funcs.periodic_commit(self._session,
                                                sequences.values()):
                del sequences[record.sequence]
                try:
                    self._realign_sequence(record.vdj, avg_len, avg_mut)

                    if record.sequence in sequences:
                        sequences[record.sequence].seq_ids += record.seq_ids
                    else:
                        sequences[record.sequence] = record
                except AlignmentException:
                    record.add_as_noresult(self._session, sample)
                except:
                    self._print('\tUnexpected error processing sequence '
                                '{}\n\t{}'.format(record.seq_ids[0],
                                                  traceback.format_exc()))

            self._print('\tAdding final sequences to database')
            for record in funcs.periodic_commit(self._session,
                                                sequences.values()):
                record.add_as_sequence(self._session, sample,
                                       meta.get('paired'))

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

def _setup_sample(session, meta):
    study, new = funcs.get_or_create(
        session, Study, name=meta.get('study_name')
    )

    if new:
        print 'Created new study "{}"'.format(study.name)
        session.commit()

    sample, new = funcs.get_or_create(
        session, Sample, name=meta.get('sample_name'),
        study=study
    )
    if new:
        sample.date = meta.get('date')
        print 'Created new sample "{}"'.format(sample.name)
        for key in ('subset', 'tissue', 'disease', 'lab',
                    'experimenter', 'ig_class', 'v_primer',
                    'j_primer'):
            setattr(sample, key, meta.get(key, require=False))
        subject, new = funcs.get_or_create(
            session, Subject, study=study,
            identifier=meta.get('subject')
        )
        sample.subject = subject
        session.commit()
    return study, sample


def get_sample_task_gen_func(session, path, trim):
    def sample_task_gen_func(task_queue):
        # Collapse identical samples
        sequence_records = {}
        parser = SeqIO.parse(
            path,
            'fasta' if path.endswith('.fasta') else 'fastq'
        )
        print 'Collapsing identical sequences'
        for record in parser:
            seq = str(record.seq)[trim:]
            if seq not in sequence_records:
                sequence_records[seq] = SequenceRecord(
                    seq, record.letter_annotations.get('phred_quality')
                )
            sequence_records[seq].seq_ids.append(record.description)
        print 'Identifying {} sequences'.format(len(sequence_records))
        for seq_record in sequence_records.values():
            task_queue.put(seq_record)

    return sample_task_gen_func


def get_align_worker_func(v_germlines, j_germlines):
    def worker_func(record, result_queue):
        try:
            vdj = VDJSequence(
                record.seq_ids[0],
                record.sequence,
                v_germlines,
                j_germlines,
                quality=record.quality
            )
            # The alignment was successful.  If the aligned sequence
            # already exists, append the seq_ids.  Otherwise add it as a
            # new unique sequence.
            record.vdj = vdj
            result_queue.put({
                'status': 'success',
                'record': record
            })
        except AlignmentException:
            result_queue.put({
                'status': 'noresult',
                'record': record
            })
        except:
            result_queue.put({
                'status': 'error',
                'message': traceback.format_exc(),
                'record': record
            })

    return worker_func


def get_collapse_func(session, sample, next_queue):
    sequence_records = {}
    def collapse_func(result):
        status = result['status']
        record = result['record']
        if status == 'success':
            if record.sequence in sequence_records:
                sequence_records[record.sequence].seq_ids.extend(
                    record.seq_ids)
            else:
                sequence_records[record.sequence] = record
                next_queue.put(result['record'])
        elif status == 'noresult':
            record.add_as_noresult(session, sample)
        elif status == 'error':
            print 'Error processing {}:\n\n{}'.format(
                record.seq_ids[0], result['message'])
    return collapse_func


def get_vties_worker_func(avg_length, avg_mut, min_similarity, max_vties):
    def vties_worker_func(record, result_queue):
        try:
            vdj = record.vdj
            vdj.align_to_germline(avg_length, avg_mut)
            if (vdj.v_match / float(vdj.v_length) < min_similarity or
                    len(vdj.v_gene) > max_vties):
                raise AlignmentException('V-match too low or too many V-ties')
            result_queue.put({
                'status': 'success',
                'record': record
            })
        except AlignmentException:
            result_queue.put({
                'status': 'noresult',
                'record': record
            })

    return vties_worker_func

def run_identify(session, args):
    mod_log.make_mod('identification', session=session, commit=True,
                     info=vars(args))
    session.close()
    # Load the germlines from files
    v_germlines = VGermlines(args.v_germlines)
    j_germlines = JGermlines(args.j_germlines, args.upstream_of_cdr3,
                             args.anchor_len, args.min_anchor_len)

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
    samples_to_process = {}
    fail = False
    for fn in sorted(metadata.keys()):
        if fn == 'all':
            continue
        meta = SampleMetadata(
            metadata[fn],
            metadata['all'] if 'all' in metadata else None
        )
        if session.query(Sample).filter(
                Sample.name == meta.get('sample_name')
                ).first() is not None:
            print 'Sample {} already exists. {}'.format(
                meta.get('sample_name'), 'Skipping.' if
                args.warn_existing else 'Cannot continue.'
            )
            fail = True
        elif meta.get('sample_name') in samples_to_process:
            print ('Sample {} exists more than once in metadata. Cannot '
                   'continue.').format(meta.get('sample_name'))
            return
        else:
            samples_to_process[meta.get('sample_name')] = {
                'path': os.path.join(args.samples_dir, fn),
                'meta': meta
            }

    if fail and not args.warn_existing:
        print ('Encountered errors.  Not running any identification.  To '
               'skip samples that are already in the database use '
               '--warn-existing.')
        return

    for sample_info in samples_to_process.values():
        # Create the study, sample, and subject if necessary
        _, sample = _setup_sample(session, sample_info['meta'])
        # Create a queue for identified sequences to be re-aligned to V-ties

        vties_queue = mp.Queue()
        # Identify sequences
        concurrent.setup_tasking(
            task_gen_func=get_sample_task_gen_func(
                session, sample_info['path'], args.trim),
            # Identify the sequences
            worker_func=get_align_worker_func(
                v_germlines=v_germlines,
                j_germlines=j_germlines,
            ),
            agg_func=get_collapse_func(session, sample, vties_queue)
        )
        vties_queue.put(None)

        # Determine average mutation rate and length
        print 'Calculating average length and mutation rate'
        realign_queue = mp.Queue()
        totals = []
        for seq in concurrent.iterate_queue(vties_queue):
            totals.append((seq.vdj.v_length, seq.vdj.mutation_fraction))
            realign_queue.put(seq)
        avg_length = sum(map(lambda r: r[0], totals)) / float(len(totals))
        avg_mut = sum(map(lambda r: r[1], totals)) / float(len(totals))
        realign_queue.put(None)

        # Realign to v-ties
        print ('Realigning {} seqs to V-ties: Avg Length={}, '
              'Avg Mut={}').format(len(totals), round(avg_length, 2),
                                   round(avg_mut, 2))
        commit_queue = mp.Queue()
        concurrent.setup_tasking_from_queue(
            realign_queue,
            worker_func=get_vties_worker_func(
                avg_length=avg_length,
                avg_mut=avg_mut,
                min_similarity=args.min_similarity / 100.0,
                max_vties=args.max_vties
            ),
            agg_func=get_collapse_func(session, sample, commit_queue)
        )
        commit_queue.put(None)

        # Finally commit the sequences
        print 'Adding final sequences to database'
        for i, record in enumerate(concurrent.iterate_queue(commit_queue)):
            record.add_as_sequence(
                session, sample, sample_info['meta'].get('paired')
            )
            if i % 100 == 0:
                print '\tCommitted {}'.format(i)
                session.commit()
        print 'Finished sample {}'.format(sample.name)
        break
