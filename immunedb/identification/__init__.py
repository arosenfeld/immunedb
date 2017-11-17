from collections import OrderedDict
import dnautils
import itertools
import traceback

from immunedb.common.models import (CDR3_OFFSET, DuplicateSequence, NoResult,
                                    Sequence)
import immunedb.util.funcs as funcs
import immunedb.util.lookups as lookups
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


def add_as_sequence(session, alignment, sample, error_action='discard'):
    try:
        seq = Sequence(
            seq_id=alignment.sequence.ids[0],
            sample_id=sample.id,

            subject_id=sample.subject.id,

            partial=alignment.partial,

            probable_indel_or_misalign=alignment.has_possible_indel,

            v_gene=funcs.format_ties(alignment.v_gene),
            j_gene=funcs.format_ties(alignment.j_gene),

            num_gaps=alignment.num_gaps,
            seq_start=alignment.seq_start,

            v_match=alignment.v_match,
            v_length=alignment.v_length,
            j_match=alignment.j_match,
            j_length=alignment.j_length,

            removed_prefix=alignment.sequence.removed_prefix_sequence,
            removed_prefix_qual=alignment.sequence.removed_prefix_quality,
            v_mutation_fraction=alignment.v_mutation_fraction,

            pre_cdr3_length=alignment.pre_cdr3_length,
            pre_cdr3_match=alignment.pre_cdr3_match,
            post_cdr3_length=alignment.post_cdr3_length,
            post_cdr3_match=alignment.post_cdr3_match,

            in_frame=alignment.in_frame,
            functional=alignment.functional,
            stop=alignment.stop,
            copy_number=len(alignment.sequence.ids),

            cdr3_nt=alignment.cdr3,
            cdr3_num_nts=len(alignment.cdr3),
            cdr3_aa=lookups.aas_from_nts(alignment.cdr3),

            sequence=str(alignment.sequence.sequence),
            quality=alignment.sequence.quality,

            locally_aligned=alignment.locally_aligned,
            insertions=alignment.insertions,
            deletions=alignment.deletions,

            germline=alignment.germline)
        session.add(seq)
        session.flush()

        # Add duplicate sequences
        try:
            session.bulk_save_objects([
                DuplicateSequence(
                    sample_id=sample.id,
                    seq_id=seq_id,
                    duplicate_seq_ai=seq.ai
                ) for seq_id in alignment.sequence.ids[1:]
            ])
        except ValueError as e:
            pass
        return seq
    except ValueError as e:
        if error_action == 'discard':
            add_as_noresult(session, alignment.sequence, sample, str(e))
            return None
        elif error_action == 'raise':
            raise e


def add_uniques(session, sample, alignments, props, aligner, realign_len=None,
                realign_mut=None):
    bucketed_seqs = OrderedDict()
    alignments = sorted(alignments, key=lambda v: v.sequence.ids[0])
    for alignment in funcs.periodic_commit(session, alignments):
        try:
            if realign_len is not None:
                aligner.align_to_germline(alignment, realign_len, realign_mut)
                if props.trim_to:
                    alignment.trim_to(props.trim_to)

            props.validate(alignment)
            bucket_key = (
                funcs.format_ties(alignment.v_gene),
                funcs.format_ties(alignment.j_gene),
                len(alignment.cdr3)
            )

            if bucket_key not in bucketed_seqs:
                bucketed_seqs[bucket_key] = {}
            bucket = bucketed_seqs[bucket_key]

            if alignment.sequence.sequence in bucket:
                bucket[alignment.sequence.sequence].sequence.ids += (
                    alignment.sequence.ids
                )
            else:
                bucket[alignment.sequence.sequence] = alignment
        except AlignmentException as e:
            add_as_noresult(session, alignment.sequence, sample, str(e))
        except Exception:
            logger.error('\tUnexpected error processing sequence '
                         '{}\n\t{}'.format(alignment.sequence.ids[0],
                                           traceback.format_exc()))

    # Collapse sequences that are the same except for Ns
    for bucket, sequences in funcs.periodic_commit(
            session, bucketed_seqs.iteritems()):
        sequences = sorted(
            sequences.values(),
            key=lambda s: (len(s.sequence.ids), s.sequence.ids[0]),
            reverse=True
        )
        while len(sequences) > 0:
            larger = sequences.pop(0)
            for i in reversed(range(len(sequences))):
                smaller = sequences[i]

                if dnautils.equal(larger.sequence.sequence,
                                  smaller.sequence.sequence):
                    larger.sequence.ids += smaller.sequence.ids
                    del sequences[i]
            add_as_sequence(session, larger, sample)
    session.commit()


def get_common_seq(seqs, cutoff=True, right=False):
    if right:
        seqs = [reversed(s) for s in seqs]
    if len(seqs) == 0:
        return seqs[0]
    gene = [
        nts[0] if all([n == nts[0] for n in nts]) else 'N'
        for nts in itertools.izip_longest(*seqs, fillvalue='N')
    ]
    if right:
        gene = reversed(gene)
    gene = ''.join(gene)
    if cutoff:
        return gene[:CDR3_OFFSET]
    return gene
