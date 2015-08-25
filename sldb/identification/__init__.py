from sldb.common.models import DuplicateSequence, NoResult, Sequence
import sldb.util.funcs as funcs
import sldb.util.lookups as lookups

class AlignmentException(Exception):
    pass


class SequenceRecord(object):
    def __init__(self, sequence, quality):
        self._orig_sequence = sequence
        self.quality = quality

        self.seq_ids = []
        self.vdj = None

    @property
    def sequence(self):
        if self.vdj is None:
            return self._orig_sequence
        return self.vdj.sequence

    def add_as_noresult(self, session, sample):
        try:
            for seq_id in self.seq_ids:
                session.add(NoResult(
                    seq_id=seq_id,
                    sample=sample,
                    sequence=self.sequence,
                    quality=funcs.ord_to_quality(self.quality)))
        except ValueError as ex:
            pass

    def add_as_sequence(self, session, sample, paired):
        try:
            indel = self.vdj.has_possible_indel
            session.add(Sequence(
                seq_id=self.vdj.id,
                sample_id=sample.id,

                paired=paired,
                partial=self.vdj.partial,

                probable_indel_or_misalign=indel,
                regions='.'.join(map(str, self.vdj.regions)),
                deletions=funcs.format_gaps(self.vdj.deletions),
                insertions=funcs.format_gaps(self.vdj.insertions),

                v_gene=funcs.format_ties(self.vdj.v_gene, 'IGHV'),
                j_gene=funcs.format_ties(self.vdj.j_gene, 'IGHJ'),

                num_gaps=self.vdj.num_gaps,
                pad_length=self.vdj.pad_length,

                v_match=self.vdj.v_match,
                v_length=self.vdj.v_length,
                j_match=self.vdj.j_match,
                j_length=self.vdj.j_length,

                removed_prefix=self.vdj.removed_prefix,
                removed_prefix_qual=funcs.ord_to_quality(
                    self.vdj.removed_prefix_qual),
                v_mutation_fraction = self.vdj.mutation_fraction,

                pre_cdr3_length=self.vdj.pre_cdr3_length,
                pre_cdr3_match=self.vdj.pre_cdr3_match,
                post_cdr3_length=self.vdj.post_cdr3_length,
                post_cdr3_match=self.vdj.post_cdr3_match,

                in_frame=self.vdj.in_frame,
                functional=self.vdj.functional,
                stop=self.vdj.stop,
                copy_number=len(self.seq_ids),

                cdr3_nt=self.vdj.cdr3,
                cdr3_num_nts=len(self.vdj.cdr3),
                cdr3_aa=lookups.aas_from_nts(self.vdj.cdr3),

                sequence=str(self.vdj.sequence),
                quality=funcs.ord_to_quality(self.vdj.quality),

                germline=self.vdj.germline))

            # Add duplicate sequences
            try:
                for seq_id in self.seq_ids:
                    if seq_id == self.vdj.id:
                        continue
                    session.add(DuplicateSequence(
                        seq_id=seq_id,
                        sample=sample,
                        duplicate_seq_id=self.vdj.id))
            except ValueError as ex:
                pass
        except ValueError as ex:
            self.add_as_noresult(session, sample)
