import dnautils
import traceback

from sldb.common.models import DuplicateSequence, NoResult, Sequence
import sldb.util.funcs as funcs
import sldb.util.lookups as lookups


class AlignmentException(Exception):
    pass


def add_as_noresult(session, vdj, sample):
    try:
        session.bulk_save_objects([
            NoResult(
                seq_id=seq_id,
                sample_id=sample.id,
                sequence=vdj.sequence,
                quality=vdj.quality
            ) for seq_id in vdj.ids
        ])
    except ValueError:
        pass


def add_as_sequence(session, vdj, sample, paired):
    try:
        seq = Sequence(
            seq_id=vdj.ids[0],
            sample_id=sample.id,

            subject_id=sample.subject.id,

            paired=paired,
            partial=vdj.partial,

            probable_indel_or_misalign=vdj.has_possible_indel,

            v_gene=funcs.format_ties(vdj.v_gene, 'IGHV'),
            j_gene=funcs.format_ties(vdj.j_gene, 'IGHJ'),

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
        except ValueError as ex:
            pass
    except ValueError:
        add_as_noresult(session, vdj, sample)


def add_uniques(session, sample, vdjs, paired, realign_len=None,
                realign_mut=None, min_similarity=None, max_vties=None):
    bucketed_seqs = {}
    for vdj in funcs.periodic_commit(session, vdjs):
        try:
            if realign_len is not None and realign_mut is not None:
                vdj.align_to_germline(realign_len, realign_mut)
                if (vdj.v_match / float(vdj.v_length) < min_similarity or
                        len(vdj.v_gene) > max_vties):
                    raise AlignmentException(
                        'V-match too low or too many V-ties'
                    )
            bucket_key = (
                funcs.format_ties(vdj.v_gene, 'IGHV'),
                funcs.format_ties(vdj.j_gene, 'IGHJ'),
                len(vdj.cdr3)
            )
            if bucket_key not in bucketed_seqs:
                bucketed_seqs[bucket_key] = {}
            bucket = bucketed_seqs[bucket_key]

            if vdj.sequence in bucket:
                bucket[vdj.sequence].ids += vdj.ids
            else:
                bucket[vdj.sequence] = vdj
        except AlignmentException:
            add_as_noresult(session, vdj, sample)
        except:
            print ('\tUnexpected error processing sequence '
                   '{}\n\t{}'.format(vdj.ids[0], traceback.format_exc()))

    # Collapse sequences that are the same except for Ns
    for sequences in funcs.periodic_commit(session, bucketed_seqs.values()):
        sequences = sorted(sequences.values(), cmp=lambda a, b:
                           cmp(len(a.ids), len(b.ids)))
        while len(sequences) > 0:
            larger = sequences.pop(0)
            for i in reversed(range(len(sequences))):
                smaller = sequences[i]

                if dnautils.equal(larger.sequence, smaller.sequence):
                    larger.ids += smaller.ids
                    del sequences[i]
            add_as_sequence(session, larger, sample, paired)
    session.commit()
