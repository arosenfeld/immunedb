import itertools

from immunedb.common.models import (CDR3_OFFSET, DuplicateSequence, NoResult,
                                    Sequence)
import immunedb.util.funcs as funcs
import immunedb.util.lookups as lookups


class AlignmentException(Exception):
    pass


def get_noresult_from_vdj(session, vdj, sample, reason):
    return [NoResult(
        seq_id=seq_id,
        sample_id=sample.id,
        sequence=vdj.orig_sequence,
        quality=vdj.orig_quality,
        reason=reason
    ) for seq_id in vdj.ids]


def get_seq_from_alignment(session, alignment, sample, strip_alleles=True):
    try:
        return [Sequence(
            seq_id=alignment.sequence.ids[0],
            sample_id=sample.id,

            subject_id=sample.subject.id,

            partial=alignment.partial,

            probable_indel_or_misalign=alignment.has_possible_indel,

            v_gene=funcs.format_ties(alignment.v_gene, strip_alleles),
            j_gene=funcs.format_ties(alignment.j_gene, strip_alleles),

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

            germline=alignment.germline)]
    except ValueError as e:
        return get_noresult_from_vdj(session, alignment.sequence, sample,
                                     str(e))


def get_duplicates_from_alignment(alignment, sample, error_action='ignore'):
    try:
        return [DuplicateSequence(
            sample_id=sample.id,
            seq_id=seq_id,
            duplicate_seq_seq_id=alignment.sequence.ids[0]
        ) for seq_id in alignment.sequence.ids[1:]]
    except ValueError as e:
        if error_action == 'ignore':
            return []
        elif error_action == 'raise':
            raise e


def add_sequences(session, alignments, sample, strip_alleles=True,
                  error_action='discard'):
    seqs_and_noresults = funcs.flatten([
        get_seq_from_alignment(session, a, sample, strip_alleles)
        for a in alignments
    ])
    succeeded = set(
        [n.seq_id for n in seqs_and_noresults if type(n) == Sequence]
    )
    funcs.bulk_add(session, seqs_and_noresults)
    session.flush()

    dups = funcs.flatten([
        get_duplicates_from_alignment(alignment, sample)
        for alignment in alignments if alignment.sequence.ids[0] in succeeded
    ])
    if dups:
        funcs.bulk_add(session, dups)


def add_noresults_for_vdj(session, vdj, sample, reason):
    try:
        funcs.bulk_add(
            session,
            get_noresult_from_vdj(session, vdj, sample, reason)
        )
    except ValueError:
        pass


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
