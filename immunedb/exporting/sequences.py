from collections import OrderedDict
from sqlalchemy.orm import joinedload

from immunedb.identification.genes import GeneName
from immunedb.exporting.writer import ExportWriter
from immunedb.common.models import (Sample, SampleMetadata, Sequence,
                                    SequenceCollapse)
from immunedb.util.log import logger
from immunedb.util.funcs import yield_limit
from immunedb.exporting.tsv_writer import StreamingTSV


mappings = {
    'airr': OrderedDict((
        ('sequence_id', 'seq_id'),
        ('sequence', 'original_sequence'),
        ('rev_comp', 'rev_comp'),
        ('productive', 'functional'),
        ('vj_in_frame', 'in_frame'),
        ('stop_codon', 'stop'),
        ('locus', lambda s: GeneName(s.v_gene).base),
        ('v_call', 'v_gene'),
        ('d_call', None),
        ('j_call', 'j_gene'),
        ('sequence_alignment', 'sequence'),
        ('germline_alignment', 'germline'),
        ('junction', 'cdr3_nt'),
        ('junction_aa', 'cdr3_aa'),
        ('v_score', 'v_match'),
        ('v_identity', lambda s: s.v_match / s.v_length),
        ('v_cigar', 'v_cigar'),
        ('d_cigar', None),
        ('j_identity', lambda s: s.j_match / s.j_length),
        ('j_score', 'j_match'),
        ('j_cigar', 'j_cigar'),
        ('v_sequence_start', 'seq_start'),
        ('duplicate_count', 'copy_number'),
        ('clone_id', 'clone_id'),

        # Extras
        ('clone_junction_nt', lambda s: s.clone.cdr3_nt if s.clone else ''),
        ('clone_junction_aa', lambda s: s.clone.cdr3_aa if s.clone else ''),
        ('sample', lambda s: s.sample.name),
        ('subject', lambda s: s.subject.identifier),
        ('duplicate_count_in_subject', lambda s:
            s.collapse.copy_number_in_subject))),


    'changeo': OrderedDict((
        ('SEQUENCE_ID', 'seq_id'),
        ('SEQUENCE_IMGT', 'sequence'),
        ('GERMLINE_IMGT_D_MASK', 'germline_d_masked'),
        ('FUNCTIONAL', 'functional'),
        ('IN_FRAME', 'in_frame'),
        ('STOP', 'stop'),
        ('V_CALL', lambda s: s.v_gene.replace(
           '|', ',' + GeneName(s.v_gene).prefix)),
        ('J_CALL', lambda s: s.j_gene.replace(
           '|', ',' + GeneName(s.j_gene).prefix)),
        ('JUNCTION_LENGTH', 'cdr3_num_nts'),
        ('JUNCTION', 'cdr3_nt'),
        ('V_IDENTITY', lambda s: round(100 * s.v_match / s.v_length, 2)),
        ('V_SCORE', 'v_match'),
        ('J_IDENTITY', lambda s: round(100 * s.j_match / s.j_length, 2)),
        ('J_SCORE', 'j_match'),
        ('DUPCOUNT', lambda s: s.copy_number),
        ('DUPCOUNT_IN_SUBJECT', lambda s: s.collapse.copy_number_in_subject),
        ('CLONE', 'clone_id'),

        # Extras
        ('CLONE_CDR3_NT', lambda s: s.clone.cdr3_nt if s.clone else ''),
        ('CLONE_CDR3_AA', lambda s: s.clone.cdr3_aa if s.clone else ''),
        ('SAMPLE', lambda s: s.sample.name),
        ('SUBJECT', lambda s: s.subject.identifier)))
}


class SequenceWriter(StreamingTSV):
    META_PREFIX = 'METADATA_'

    def __init__(self, format_name, metadata_fields):
        self.mapping = mappings[format_name]
        fields = list(self.mapping.keys()) + list(sorted([
            self.META_PREFIX + m for m in metadata_fields]))
        super(SequenceWriter, self).__init__(fields)

    def writeseq(self, seq):
        return super(SequenceWriter, self).writerow(self.format_seq(seq))

    def format_seq(self, seq):
        fields = {
            out_field: self._get_val(seq, model_field)
            for out_field, model_field in self.mapping.items()
        }

        fields.update({
            self.META_PREFIX + k: v
            for k, v in seq.sample.metadata_dict.items()
        })
        return fields

    def _get_val(self, model, field):
        if type(field) == str:
            val = getattr(model, field)
            if type(val) == bool:
                return 'T' if val else 'F'
            return val
        if callable(field):
            return field(model)
        return ''


def get_sequences(session, sample, fmt, clones_only, min_subject_copies):
    meta_keys = set([m.key for m in session.query(SampleMetadata.key)])

    seqs = session.query(Sequence).filter(
        Sequence.sample_id == sample.id
    ).join(
        SequenceCollapse
    ).options(
        joinedload(Sequence.clone),
        joinedload(Sequence.collapse),
        joinedload(Sequence.sample),
        joinedload(Sequence.subject),
    )

    if clones_only:
        seqs = seqs.filter(~Sequence.clone_id.is_(None))
    if min_subject_copies:
        seqs = seqs.filter(
            SequenceCollapse.copy_number_in_subject >=
            min_subject_copies
        )

    writer = SequenceWriter(fmt, meta_keys)

    yield writer.writeheader()
    for seq in yield_limit(seqs, Sequence.ai):
        yield writer.writeseq(seq)


def write_sequences(session, sample_ids=None, out_format='changeo',
                    clones_only=False, min_subject_copies=None, zipped=False,
                    **kwargs):
    samples = session.query(Sample)
    if sample_ids:
        samples = samples.filter(Sample.id.in_(sample_ids))
    with ExportWriter(zipped=zipped) as fh:
        for sample in samples:
            logger.info('Exporting sample {}'.format(sample.name))
            fh.set_filename('{}.{}.tsv'.format(sample.name, out_format))
            fh.write(
                get_sequences(session, sample, out_format, clones_only,
                              min_subject_copies)
            )
        return fh.get_zip_value()
