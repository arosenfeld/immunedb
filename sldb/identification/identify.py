import distance
import json
import os
import re
import time

from Bio import SeqIO
from Bio.Seq import Seq

from sldb.common.config import allowed_read_types
from sldb.identification.vdj_sequence import VDJSequence
from sldb.identification.v_genes import VGermlines
from sldb.common.models import *
import sldb.util.funcs as funcs
import sldb.util.lookups as lookups


class SampleMetadata(object):
    def __init__(self, global_config, specific_config):
        self._global = global_config
        self._specific = specific_config

    def get(self, key, require=True):
        if key in self._specific:
            return self._specific[key]
        if key in self._global:
            return self._global[key]
        if require:
            raise Exception(('Could not find metadata for key {}'.format(key)))



def add_to_db(session, alignment, sample, vdj):
    existing = session.query(Sequence).filter(
        Sequence.sequence == vdj.sequence,
        Sequence.sample == sample).first()
    if existing is not None:
        existing.copy_number += 1
        session.add(DuplicateSequence(duplicate_seq_id=existing.seq_id,
                                      sample_id=sample.id,
                                      seq_id=vdj.id))
    else:
        session.add(Sequence(
            seq_id=vdj.id,
            sample=sample,

            alignment=alignment,
            probable_indel_or_misalign=vdj.has_possible_indel,

            v_gene=funcs.format_ties(vdj.v_gene, 'IGHV'),
            j_gene=funcs.format_ties(vdj.j_gene, 'IGHJ'),

            num_gaps=vdj.num_gaps,
            pad_length=vdj.pad_length,

            v_match=vdj.v_match,
            v_length=vdj.v_length,
            j_match=vdj.j_match,
            j_length=vdj.j_length,

            pre_cdr3_length=vdj.pre_cdr3_length,
            pre_cdr3_match=vdj.pre_cdr3_match,
            post_cdr3_length=vdj.post_cdr3_length,
            post_cdr3_match=vdj.post_cdr3_match,

            in_frame=vdj.in_frame,
            functional=vdj.functional,
            stop=vdj.stop,
            copy_number=vdj.copy_number,

            junction_nt=vdj.cdr3,
            junction_aa=lookups.aas_from_nts(vdj.cdr3, ''),
            gap_method='IMGT',

            sequence=str(vdj.sequence),
            sequence_replaced=vdj.sequence_filled,

            germline=vdj.germline))


def _identify_reads(session, path, fn, meta, v_germlines, limit_alignments):
    print 'Starting {}'.format(fn)
    study, new = funcs.get_or_create(
        session, Study, name=meta.get('study_name'))

    read_type = meta.get('read_type')
    if read_type not in allowed_read_types:
        raise Exception(('Invalid read type {}.  Must be'
                         'one of {}').format(
                        read_type, ','.join(allowed_read_types)))
    if read_type not in limit_alignments:
        print ('Skipping since read type is {} and identification '
               'limited to {}').format(read_type, ','.join(limit_alignments))
        return

    if new:
        print 'Created new study "{}"'.format(study.name)
        session.commit()
    else:
        print 'Study "{}" already exists in MASTER'.format(study.name)

    sample, new = funcs.get_or_create(
        session, Sample, name=meta.get('sample_name', require=False) or
        fn.split('.', 1)[0], study=study)
    if new:
        sample.date = meta.get('date')
        print '\tCreated new sample "{}" in MASTER'.format(sample.name)
        for key in ('subset', 'tissue', 'disease', 'lab',
                    'experimenter'):
            setattr(sample, key, meta.get(key, require=False))
        subject, new = funcs.get_or_create(
            session, Subject, study=study, identifier=meta.get('subject'))
        sample.subject = subject
        session.commit()
    else:
        # TODO: Verify the data for the existing sample matches
        exists = session.query(Sequence).filter(
            Sequence.sample == sample,
            Sequence.alignment == read_type).first()
        if exists is not None:
            print ('\tSample "{}" for study already exists in DATA.  '
                   'Skipping.').format(sample.name)
            return

    # Sums of statistics for retroactive v-tie calculation
    lengths_sum = 0
    mutations_sum = 0
    # All VDJs assigned for this sample, keyed by raw sequence
    vdjs = {}
    # All probable duplicates based on keys of vdjs dictionary
    dups = {}
    # Number of noresult sequences
    no_result = 0

    start = time.time()
    print '\tIdentifying V and J, committing No Results'
    for i, record in enumerate(SeqIO.parse(os.path.join(path, fn), 'fasta')):
        if i > 0 and i % 1000 == 0:
            print '\t\tCommitted {}'.format(i)
            session.commit()
        # Key the vdjs dictionary by the unmodified sequence
        key = str(record.seq)
        if key in vdjs:
            # If this exact sequence, without padding or gaps, has been
            # assigned a V and J, bump the copy number of that
            # VDJSequence instance and add this seq_id as a duplicate.
            vdj = vdjs[key]
            vdj.copy_number += 1
            if vdj.id not in dups:
                dups[vdj.id] = []
            dups[vdj.id].append(DuplicateSequence(
                duplicate_seq_id=vdjs[key].id,
                sample_id=sample.id,
                seq_id=record.description))
        else:
            # This is the first instance of this exact sequence, so align it
            # and identify it's V and J
            vdj = VDJSequence(record.description,
                              record.seq,
                              read_type == 'R1+R2',
                              v_germlines)
            if vdj.v_gene is not None and vdj.j_gene is not None:
                # If the V and J are found, add it to the vdjs dictionary to
                # prevent future exact copies from being aligned
                lengths_sum += vdj.v_length
                mutations_sum += vdj.mutation_fraction
                vdjs[key] = vdj
            else:
                # The V or J could not be found, so add it as a noresult
                session.add(NoResult(sample=sample,
                                     seq_id=vdj.id,
                                     sequence=str(vdj.sequence)))
                no_result += 1

    session.commit()
    print '\t\tCumulative time: {} sec'.format(round(time.time() - start))

    if len(vdjs) == 0:
        print '\t\tNo sequences identified'
        return

    print '\tCalculating V-ties'
    avg_len = lengths_sum / float(len(vdjs))
    avg_mut = mutations_sum / float(len(vdjs))

    for i, (_, vdj) in enumerate(vdjs.iteritems()):
        if i > 0 and i % 1000 == 0:
            print '\t\tCommitted {}'.format(i)
            session.commit()
        # Align the sequence to a germline based on v_ties
        vdj.align_to_germline(avg_len, avg_mut)
        if vdj.v_gene is not None and vdj.j_gene is not None:
            # Add the sequence to the database
            add_to_db(session, read_type, sample, vdj)
        else:
            # This is a rare condition, but some sequence after aligning to
            # V-ties the CDR3 becomes non-existent, and it is thrown out
            no_result += 1
            session.add(NoResult(sample=sample,
                                 seq_id=vdj.id,
                                 sequence=str(vdj.sequence)))
            if vdj.id in dups:
                # It's duplicate sequences must be added as noresults also
                for dup in dups[vdj.id]:
                    # Restore the original sequence by removing padding and
                    # gaps
                    session.add(NoResult(
                        sample=sample,
                        seq_id=dup.seq_id,
                        sequence=vdj.sequence.replace('-', '').strip('N')))

                del dups[vdj.id]
    session.commit()

    # Add the true duplicates to the database
    print '\tAdding duplicates'
    for i, dup_seqs in enumerate(dups.values()):
        if i > 0 and i % 1000 == 0:
            print '\t\tCommitted {}'.format(i)
            session.commit()
        session.add_all(dup_seqs)
    session.commit()

    print '\t\tCumulative time: {} sec'.format(round(time.time() - start))

    print '\tlen={}'.format(avg_len)
    print '\tmut={}'.format(avg_mut)
    print '\tnoresults={}'.format(no_result)


def run_identify(session, args):
    v_germlines = VGermlines(args.v_germlines)

    for base_dir in args.base_dirs:
        print 'Descending into {}'.format(base_dir)
        meta_fn = '{}/metadata.json'.format(base_dir)
        if not os.path.isfile(meta_fn):
            print '\tNo metadata file found for this set of samples.'
            return

        with open(meta_fn) as fh:
            metadata = json.load(fh)
            for fn in os.listdir(base_dir):
                if fn == 'metadata.json' or fn not in metadata:
                    continue
                _identify_reads(
                    session,
                    base_dir,
                    fn,
                    SampleMetadata(metadata[fn], metadata['all']),
                    v_germlines,
                    args.limit_alignments)
