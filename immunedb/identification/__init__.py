from collections import OrderedDict
import dnautils
import itertools
import re
import traceback

from immunedb.common.models import (CDR3_OFFSET, DuplicateSequence, NoResult,
                                    Sequence)
import immunedb.util.funcs as funcs
import immunedb.util.lookups as lookups
from immunedb.util.hyper import hypergeom
from immunedb.util.log import logger


class AlignmentException(Exception):
    pass


def add_as_noresult(session, vdj, sample, reason):
    try:
        session.bulk_save_objects([
            NoResult(
                seq_id=seq_id,
                sample_id=sample.id,
                sequence=vdj.orig_sequence,
                quality=vdj.orig_quality,
                reason=reason
            ) for seq_id in vdj.ids
        ])
    except ValueError:
        pass


def add_as_sequence(session, vdj, sample):
    try:
        seq = Sequence(
            seq_id=vdj.ids[0],
            sample_id=sample.id,

            subject_id=sample.subject.id,

            partial=vdj.partial,

            probable_indel_or_misalign=vdj.has_possible_indel,

            v_gene=funcs.format_ties(vdj.v_gene),
            j_gene=funcs.format_ties(vdj.j_gene),

            num_gaps=vdj.num_gaps,
            pad_length=vdj.pad_length,

            v_match=vdj.v_match,
            v_length=vdj.v_length,
            j_match=vdj.j_match,
            j_length=vdj.j_length,

            removed_prefix=vdj.removed_prefix,
            removed_prefix_qual=vdj.removed_prefix_qual,
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
            quality=vdj.quality,

            germline=vdj.germline)
        session.add(seq)
        session.flush()

        # Add duplicate sequences
        try:
            session.bulk_save_objects([
                DuplicateSequence(
                    sample_id=sample.id,
                    seq_id=seq_id,
                    duplicate_seq_ai=seq.ai
                ) for seq_id in vdj.ids[1:]
            ])
        except ValueError:
            pass
    except ValueError as e:
        add_as_noresult(session, vdj, sample, str(e))


def add_uniques(session, sample, vdjs, props, realign_len=None,
                realign_mut=None):
    bucketed_seqs = OrderedDict()
    vdjs = sorted(vdjs, key=lambda v: v.ids[0])
    for vdj in funcs.periodic_commit(session, vdjs):
        try:
            if realign_len is not None:
                vdj.align_to_germline(realign_len, realign_mut, props.trim_to)
            if not props.valid_min_similarity(vdj):
                raise AlignmentException(
                    'V-identity too low {} < {}'.format(
                        vdj.v_match / float(vdj.v_length),
                        props.min_similarity))
            if not props.valid_v_ties(vdj):
                raise AlignmentException('Too many V-ties {} > {}'.format(
                    len(vdj.v_gene), props.max_v_ties))
            if not props.valid_padding(vdj):
                raise AlignmentException('Too much padding {} (max {})'.format(
                    vdj.pad_length, props.max_padding))
            if not props.valid_families(vdj):
                raise AlignmentException('Cross-family V-call')
            bucket_key = (
                funcs.format_ties(vdj.v_gene),
                funcs.format_ties(vdj.j_gene),
                len(vdj.cdr3)
            )

            if bucket_key not in bucketed_seqs:
                bucketed_seqs[bucket_key] = {}
            bucket = bucketed_seqs[bucket_key]

            if vdj.sequence in bucket:
                bucket[vdj.sequence].ids += vdj.ids
            else:
                bucket[vdj.sequence] = vdj
        except AlignmentException as e:
            add_as_noresult(session, vdj, sample, str(e))
        except:
            logger.error('\tUnexpected error processing sequence '
                         '{}\n\t{}'.format(vdj.ids[0], traceback.format_exc()))

    # Collapse sequences that are the same except for Ns
    for bucket, sequences in funcs.periodic_commit(
            session, bucketed_seqs.iteritems()):
        sequences = sorted(
            sequences.values(),
            key=lambda s: (len(s.ids), s.ids[0]),
            reverse=True
        )
        while len(sequences) > 0:
            larger = sequences.pop(0)
            for i in reversed(range(len(sequences))):
                smaller = sequences[i]

                if dnautils.equal(larger.sequence, smaller.sequence):
                    larger.ids += smaller.ids
                    del sequences[i]
            add_as_sequence(session, larger, sample)
    session.commit()


def get_common_seq(seqs, cutoff=True):
    if len(seqs) == 0:
        return seqs[0]
    v_gene = []
    for nts in itertools.izip_longest(*seqs, fillvalue='N'):
        v_gene.append(nts[0] if all(map(lambda n: n == nts[0], nts)) else 'N')
    v_gene = ''.join(v_gene)
    if cutoff:
        return v_gene[:CDR3_OFFSET]
    return v_gene
