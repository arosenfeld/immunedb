import distance
import json
import os
import re

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


def _create_mapping(session, identity_seq_id, alignment, sample, vdj):
    if vdj.possible_indel:
        seq = vdj.sequence.lstrip('N')
        germ = vdj.germline[len(vdj.sequence) - len(seq):]
        lev_dist = distance.levenshtein(germ, seq)
    else:
        lev_dist = None

    return SequenceMapping(
        identity_seq_id=identity_seq_id,
        sample=sample,
        seq_id=vdj.id,

        alignment=alignment,
        levenshtein_dist=lev_dist,

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
        stop=vdj.stop,
        copy_number=1,
        functional=vdj.functional,

        sequence=str(vdj.sequence))


def _format_ties(ties, name):
    if ties is None:
        return None
    ties = map(lambda e: e.replace(name, ''), ties)
    return '{}{}'.format(name, '|'.join(sorted(ties)))


def _add_to_db(session, alignment, sample, vdj):
    # Check if a sequence with the exact same filled sequence exists
    m = session.query(Sequence).filter(
        Sequence.sequence_replaced == str(vdj.sequence_filled)).first()
    if m is None:
        # If not, add the sequence, and make the mapping
        m = Sequence(seq_id=vdj.id,
                     v_call=_format_ties(vdj.v_gene, 'IGHV'),
                     j_call=_format_ties(vdj.j_gene, 'IGHJ'),
                     junction_num_nts=len(vdj.cdr3),
                     junction_nt=str(vdj.cdr3),
                     junction_aa=lookups.aas_from_nts(vdj.cdr3, ''),
                     gap_method='IMGT',
                     sequence_replaced=str(vdj.sequence_filled),
                     germline=vdj.germline)
        session.add(m)
        session.add(_create_mapping(session, m.seq_id, alignment, sample, vdj))
    else:
        # If there is an identical sequence, check if its appeared in this
        # sample.
        existing = session.query(SequenceMapping).filter(
            SequenceMapping.unique_id == \
                funcs.hash(m.seq_id, sample.id, vdj.sequence)).first()

        if existing is not None:
            # If so, bump the copy number and insert the duplicate sequence
            existing.copy_number += 1
            session.add(DuplicateSequence(duplicate_seq_id=existing.seq_id,
                                          sample_id=sample.id,
                                          seq_id=vdj.id))
        else:
            # If not, add a new mapping from the existing sequence
            session.add(_create_mapping(session, m.seq_id, alignment, sample,
                                        vdj))


def _identify_reads(session, path, fn, meta, v_germlines, full_only):
    print 'Starting {}'.format(fn)
    study, new = funcs.get_or_create(
        session, Study, name=meta.get('study_name'))

    read_type = meta.get('read_type')
    if read_type not in allowed_read_types:
        raise Exception(('Invalid read type {}.  Must be'
                         'one of {}').format(
                            read_type, ','.join(allowed_read_types)))
    if read_type == 'R1+R2' and full_only:
        print ('Skpping {} since it contains partial reads and read_type '
               'is "R1+R2"'.format(sample.name))
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
        print '\tCreated new sample "{}" in MASTER'.format(sample.name)
        for key in ('date', 'subset', 'tissue', 'disease', 'lab',
                    'experimenter'):
            setattr(sample, key, meta.get(key, require=False))
        subject, new = funcs.get_or_create(
            session, Subject, study=study, identifier=meta.get('subject'))
        sample.subject = subject
        session.commit()
    else:
        # TODO: Verify the data for the existing sample matches
        exists = session.query(SequenceMapping).filter(
            SequenceMapping.sample == sample,
            SequenceMapping.alignment == read_type).first()
        if exists is not None:
            print ('\tSample "{}" for study already exists in DATA.  '
                   'Skipping.').format(sample.name)
            return

    lengths_sum = 0
    mutations_sum = 0
    vdjs = []
    no_result = 0

    for i, record in enumerate(SeqIO.parse(os.path.join(path, fn), 'fasta')):
        if i > 0 and i % 1000 == 0:
            print '\tCommitted {}'.format(i)
            session.commit()

        vdj = VDJSequence(record.description, 
                          record.seq,
                          read_type == 'R1+R2',
                          v_germlines)
        if vdj.v_gene is not None and vdj.j_gene is not None:
            lengths_sum += vdj.v_length
            mutations_sum += vdj.mutation_fraction
            vdjs.append(vdj)
        else:
            session.add(NoResult(sample=sample,
                                 seq_id=vdj.id,
                                 sequence=str(vdj.sequence)))
            no_result += 1
    session.commit()
    cnt = i

    if len(vdjs) == 0:
        print '\tNo sequences identified'
        return

    avg_len = lengths_sum / float(len(vdjs))
    avg_mut = mutations_sum / float(len(vdjs))

    for vdj in vdjs:
        vdj.align_to_germline(avg_len, avg_mut)
        _add_to_db(session, read_type, sample, vdj)

    print '\tlen={}'.format(avg_len)
    print '\tmut={}'.format(avg_mut)
    print '\tnoresults={}'.format(no_result)
    session.commit()


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
                if fn == 'metadata.json':
                    continue
                _identify_reads(
                    session,
                    base_dir,
                    fn,
                    SampleMetadata(metadata[fn], metadata['all']),
                    v_germlines,
                    args.only_full)
