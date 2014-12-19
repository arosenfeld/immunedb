import distance
import json
import os
import re
from Bio import SeqIO

from sldb.identification.vdj_sequence import VDJSequence
from sldb.identification.v_ties import get_v_ties
from sldb.common.models import *
import sldb.util.funcs as funcs
import sldb.util.lookups as lookups


def _get_from_meta(meta, sample_name, key, require):
    if sample_name in meta and key in meta[sample_name]:
        return meta[sample_name][key]
    if key in meta['all']:
        return meta['all'][key]
    if require:
        raise Exception(('Unknown key {} '
                        'for sample {}').format(key, sample_name))
    return None


def _create_mapping(session, identity_seq_id, alignment, sample, vdj):
    if vdj.probable_deletion or vdj.v_match / float(vdj.v_length) < .6:
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


def _add_to_db(session, alignment, sample, vdj):
    assert alignment in ['R1', 'R2', 'pRESTO']
    # Check if a sequence with the exact same filled sequence exists
    m = session.query(Sequence).filter(
        Sequence.sequence_replaced == str(vdj.sequence_filled)).first()
    if m is None:
        # If not, add the sequence, and make the mapping
        m = Sequence(seq_id=vdj.id,
                     v_call='|'.join(vdj.v_gene),
                     j_call=vdj.j_gene,
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
            if session.query(DuplicateSequence).filter(
                    DuplicateSequence.identity_seq_id== m.seq_id,
                    DuplicateSequence.seq_id == vdj.id).first() == None:
                session.add(DuplicateSequence(identity_seq_id=m.seq_id,
                                              seq_id=vdj.id))
        else:
            # If not, add a new mapping from the existing sequence
            session.add(_create_mapping(session, m.seq_id, alignment, sample,
                                        vdj))


def _identify_reads(session, meta, fn, name, alignment):
    print 'Starting {}'.format(fn)
    study, new = funcs.get_or_create(session, Study, name=meta['all']['study'])
    if new:
        print 'Created new study "{}"'.format(study.name)
        session.commit()
    else:
        print 'Study "{}" already exists in MASTER'.format(study.name)

    sample, new = funcs.get_or_create(session, Sample, name=name, study=study)
    if new:
        print '\tCreated new sample "{}" in MASTER'.format(sample.name)
        for key in ['date', 'subset', 'tissue', 'disease', 'lab',
                    'experimenter']:
            setattr(sample, key, _get_from_meta(meta, name, key, False))
        subject, new = funcs.get_or_create(
            session, Subject, study=study, identifier=_get_from_meta(
                meta, name, 'subject', True))
        sample.subject = subject
        session.commit()
    else:
        # TODO: Verify the data for the existing sample matches
        exists = session.query(SequenceMapping).filter(
            SequenceMapping.sample == sample,
            SequenceMapping.alignment == alignment).first()
        if exists is not None:
            print ('\tSample "{}" for study already exists in DATA.  '
                   'Skipping.').format(sample.name)
            return

    lengths_sum = 0
    mutations_sum = 0
    vdjs = []
    no_result = 0

    for i, record in enumerate(SeqIO.parse(fn, 'fasta')):
        if i > 0 and i % 1000 == 0:
            print '\tCommitted {}'.format(i)
            session.commit()

        vdj = VDJSequence(record.description, record.seq, alignment == 'pRESTO')
        if vdj.v_gene is not None and vdj.j_gene is not None:
            lengths_sum += vdj.v_length
            mutations_sum += vdj.mutation_fraction
            vdjs.append(vdj)
            _add_to_db(session, alignment, sample, vdj)
        else:
            session.add(NoResult(sample=sample, seq_id=vdj.id))
            no_result += 1
    cnt = i

    if len(vdjs) == 0:
        print '\tNo sequences identified'
        return

    avg_len = lengths_sum / float(len(vdjs))
    avg_mutation_frac = mutations_sum / float(len(vdjs))
    v_ties = get_v_ties(avg_len, avg_mutation_frac)

    v_tie_cnt = 0
    for vdj in vdjs:
        new_vs = set([])
        for v in vdj.v_gene:
            if v in v_ties:
                ties = v_ties[v].split('|')
                new_vs = new_vs.union(set(ties))
            else:
                new_vs.add(v)
        old_vs = vdj.v_gene[:]
        vdj.v_gene = list(new_vs)
        if len(vdj.v_gene) > 1:
            v_tie_cnt += 1
    print '\ttotal_seqs={}'.format(cnt)
    print '\tvties={}'.format(v_tie_cnt)
    print '\tlen={}'.format(avg_len)
    print '\tmut={}'.format(avg_mutation_frac)
    print '\tnoresults={}'.format(no_result)


def run_identify(session, args):
    for base_dir in args.base_dirs: 
        print 'Parsing {}'.format(base_dir)
        meta_fn = '{}/metadata.json'.format(base_dir)
        if not os.path.isfile(meta_fn):
            print 'No metadata file found for this set of samples.'
            return
        with open(meta_fn) as fh:
            metadata = json.load(fh)

            names = set([])
            for fn in os.listdir('{}/processed'.format(base_dir)):
                name = fn.split('.')[0].rsplit('_', 2)[0]
                names.add(name)

            for name in names:
                base = '{}/presto/{}'.format(base_dir, name)
                r1 = '{}_R1_001.sync_assemble-fail.fasta'.format(base)
                r2 = '{}_R2_001.sync_assemble-fail.fasta'.format(base)
                join = '{}_R1_001.sync_assemble-pass.fasta'.format(base)
                if not os.path.isfile('{}.log'.format(base)):
                    print 'Skipping {} since no presto log exists.'.format(name)
                    continue
                _identify_reads(session, metadata,
                                join, base.rsplit('/', 1)[1], 'pRESTO')
                _identify_reads(session, metadata,
                                r1, base.rsplit('/', 1)[1], 'R1')
