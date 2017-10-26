from collections import OrderedDict

from immunedb.common.models import Sequence, SequenceCollapse
from immunedb.exporting.export import Exporter


def _if_clone(f):
    def _wrap(seq):
        if seq.clone_id is not None:
            return f(seq)
        return ''
    return _wrap


class SequenceExport(Exporter):
    """A class to handle exporting sequences in various formats.

    :param Session session: A handle to a database session
    :param str eformat: The export format to use.  Currently "csv", "orig", and
        "clip" for tab-delimited, FASTA, FASTA with filled in germlines, and
        FASTA in CLIP format respectively
    :param str rtype: The type of record to filter the query on.  Currently
        either "sample" or "clone"
    :param int-array rids: A list of IDs of ``rtype`` to export
    :param str-array selected_fields: A list of fields to export
    :param bool noresults: If sequences which could not be assigned a V or J
        should be included in the output

    """

    allowed_fields = OrderedDict({
        'seq_id': lambda seq: seq.seq_id,
        'subject_id': lambda seq: seq.sample.subject.id,
        'subject_identifier': lambda seq: seq.sample.subject.identifier,
        'subset': lambda seq: seq.sample.subset,
        'ig_class': lambda seq: seq.sample.ig_class,
        'tissue': lambda seq: seq.sample.tissue,
        'disease': lambda seq: seq.sample.disease,
        'lab': lambda seq: seq.sample.lab,
        'experimenter': lambda seq: seq.sample.experimenter,
        'v_primer': lambda seq: seq.sample.v_primer,
        'j_primer': lambda seq: seq.sample.j_primer,
        'date': lambda seq: seq.sample.date,

        'sample_id': lambda seq: seq.sample_id,
        'sample_name': lambda seq: seq.sample.name,
        'study_id': lambda seq: seq.sample.study.id,
        'study_name': lambda seq: seq.sample.study.name,

        'partial': lambda seq: seq.partial,
        'probable_indel_or_misalign':
            lambda seq: seq.probable_indel_or_misalign,
        'num_gaps': lambda seq: seq.num_gaps,
        'seq_start': lambda seq: seq.seq_start,
        'v_match': lambda seq: seq.v_match,
        'v_length': lambda seq: seq.v_length,
        'j_match': lambda seq: seq.j_match,
        'j_length': lambda seq: seq.j_length,
        'pre_cdr3_match': lambda seq: seq.pre_cdr3_match,
        'pre_cdr3_length': lambda seq: seq.pre_cdr3_length,
        'post_cdr3_match': lambda seq: seq.post_cdr3_match,
        'post_cdr3_length': lambda seq: seq.post_cdr3_length,
        'in_frame': lambda seq: seq.in_frame,
        'functional': lambda seq: seq.functional,
        'stop': lambda seq: seq.stop,
        'copy_number': lambda seq: seq.copy_number,
        'sequence': lambda seq: seq.sequence,
        'quality': lambda seq: seq.quality,
        'germline': lambda seq: seq.germline,
        'v_gene': lambda seq: seq.v_gene,
        'j_gene': lambda seq: seq.j_gene,
        'cdr3_nt': lambda seq: seq.cdr3_nt,
        'cdr3_aa': lambda seq: seq.cdr3_aa,
        'cdr3_num_nts': lambda seq: seq.cdr3_num_nts,

        'clone_id': _if_clone(lambda seq: seq.clone_id),
        'clone_cdr3_nt': _if_clone(lambda seq: seq.clone.cdr3_nt),
        'clone_cdr3_aa': _if_clone(lambda seq: seq.clone.cdr3_aa),
        'clone_json_tree': _if_clone(lambda seq: seq.clone.tree),
        'clone_parent_id': _if_clone(lambda seq: seq.clone.parent_id),

        'copy_number_in_subject':
            lambda seq: seq.collapse.copy_number_in_subject,
        'collapse_to_subject_sample_id': lambda seq:
            seq.collapse.collapse_to_subject_sample_id,
        'collapse_to_subject_seq_id': lambda seq:
            seq.collapse.collapse_to_subject_seq_id,
        'instances_in_subject': lambda seq:
            seq.collapse.instances_in_subject,
    })

    def __init__(self, session, writer, rtype, rids, selected_fields,
                 subject_uniques, only_with_clones):
        super(SequenceExport, self).__init__(
            session, rtype, rids, SequenceExport.allowed_fields,
            selected_fields)
        self.writer = writer
        self.subject_uniques = subject_uniques
        self.only_with_clones = only_with_clones

    def get_data(self):
        """Gets the output data for this export instance.  This could be export
        to a file, over a socket, etc. as it's simply a string or could be
        treated as a byte string.

        :returns: Lines of output for the selected sequences
        :rtype: str

        """

        # Get all the sequences matching the request
        self.writer.set_selected_fields(self.selected_fields)
        queries = [self.session.query(Sequence).filter(
            getattr(
                Sequence, '{}_id'.format(self.rtype)
            ) == rid
        ) for rid in sorted(self.rids)]
        for query in queries:
            query = query.join(SequenceCollapse)
            if self.subject_uniques:
                query = query.filter(
                    SequenceCollapse.copy_number_in_subject >= 1
                )
            if self.only_with_clones:
                query = query.filter(~Sequence.clone_id.is_(None))

            query = self.writer.preprocess_query(query)
            self.selected_fields = self.writer.get_required_fields(
                self.selected_fields)

            for seq in query:
                # Get the selected data for the sequence
                data = self.get_selected_data(seq)
                yield self.writer.format_sequence(data)
