from sldb.common.models import Sequence, DuplicateSequence, NoResult
from sldb.exporting.export import Exporter
from sldb.util.nested_writer import NestedCSVWriter


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
    :param bool min_copy: Minimum copy number of sequences to include
    :param bool duplicates: If duplicate sequences should be included in the
        output
    :param bool noresults: If sequences which could not be assigned a V or J
        should be included in the output

    """
    _allowed_fields = [
        'seq_id',
        ('duplicate_of_seq_id', lambda seq: None),
        ('subject_id', lambda seq: seq.sample.subject.id),
        ('subject_identifier', lambda seq: seq.sample.subject.identifier),
        ('subset', lambda seq: seq.sample.subset),
        ('tissue', lambda seq: seq.sample.tissue),
        ('disease', lambda seq: seq.sample.disease),
        ('lab', lambda seq: seq.sample.lab),
        ('experimenter', lambda seq: seq.sample.experimenter),
        ('date', lambda seq: seq.sample.date),

        'sample_id',
        ('sample_name', lambda seq: seq.sample.name),
        ('study_id', lambda seq: seq.sample.study.id),
        ('study_name', lambda seq: seq.sample.study.name),

        'alignment', 'partial', 'probable_indel_or_misalign',

        'num_gaps', 'pad_length',

        'v_match', 'v_length', 'j_match', 'j_length',

        'pre_cdr3_match', 'pre_cdr3_length', 'post_cdr3_match',
        'post_cdr3_length',

        'in_frame', 'functional', 'stop', 'copy_number',

        'sequence', 'quality', 'germline',

        'v_gene', 'j_gene', 'cdr3_nt', 'cdr3_aa', 'cdr3_num_nts', 'gap_method',

        'clone_id',
        ('clone_cdr3_nt', lambda seq: seq.clone.cdr3_nt),
        ('clone_cdr3_aa', lambda seq: seq.clone.cdr3_aa),
        ('clone_json_tree', lambda seq: seq.clone.tree),

        'copy_number_in_sample', 'collapse_to_sample_seq_id',

        'copy_number_in_subject', 'collapse_to_subject_seq_id',
        'collapse_to_subject_sample_id',
    ]

    def __init__(self, session, writer, rtype, rids, selected_fields,
                 min_copy, duplicates, noresults, level, only_with_clones):
        super(SequenceExport, self).__init__(
            session, rtype, rids, SequenceExport._allowed_fields,
            selected_fields)
        self.writer = writer
        self.min_copy = min_copy
        self.duplicates = duplicates
        self.noresults = noresults
        self.level = level
        self.only_with_clones = only_with_clones

    def _fasta_entry(self, seq_id, info, sequence):
        """Gets the entry for a sequence in FASTA format with ``info`` separated
        in the header by pipe symbols.

        :param str seq_id: The sequence identifier
        :param dict info: A dictionary of info to include in the header
        :param str sequence: The sequence for the entry

        :returns: The FASTA entry like

            >seq_id|hdr_key1=hdr_val1|hdr_key2=hdr_val2|...
            ATCGATCG...

        :rtype: str

        """
        return '>{}{}{}\n{}\n'.format(
            seq_id,
            '|' if len(info) > 0 else '',
            '|'.join(map(
                lambda (k, v): '{}={}'.format(k, v), info.iteritems())),
            sequence)

    def get_data(self):
        """Gets the output data for this export instance.  This could be export
        to a file, over a socket, etc. as it's simply a string or could be
        treated as a byte string.

        :returns: Lines of output for the selected sequences
        :rtype: str

        """
        # Get all the sequences matching the request
        seqs = self.session.query(Sequence).filter(
            getattr(
                Sequence, '{}_id'.format(self.rtype)
            ).in_(self.rids)
        )
        if self.level == 'uncollapsed':
            seqs = seqs.filter(Sequence.copy_number >= self.min_copy)
        else:
            seqs = seqs.filter(
                getattr(Sequence, 'copy_number_in_{}'.format(self.level)) >=
                self.min_copy
            )

        if self.only_with_clones:
            seqs = seqs.filter(~Sequence.clone_id.is_(None))

        # If it's a CLIP file, order by clone_id to minimize
        # repetition of germline entries
        seqs = self.writer.preprocess_query(seqs)
        headers = []
        for field in self.export_fields:
            n, f = self._name_and_field(field)
            if n in self.selected_fields:
                headers.append(n)
        headers = self.writer.preprocess_headers(headers)

        # For CLIP files to check if the germline needs to be output
        last_cid = ''
        # TODO: Change this to .yield_per?
        for seq in seqs:
            # Get the selected data for the sequence
            data = self.get_selected_data(seq)
            yield self.writer.format_sequence(data)

            # Regardless of output type, output duplicates if desired
            if self.duplicates and seq.copy_number > 1:
                for dup in self.session.query(DuplicateSequence).filter(
                        DuplicateSequence.duplicate_seq == seq):
                    # Duplicates have the same data as their original sequence,
                    # except the seq_id, copy_number is set to 0, and the
                    # duplicate_of_seq_id is set to the original sequence's
                    data = self.get_selected_data(
                        seq,
                        seq_id=dup.seq_id,
                        copy_number=0,
                        duplicate_of_seq_id=seq.seq_id)

                    yield self.writer.format_sequence(data, pseudo_seq=True)

        # Regardless of the output type, output noresults if desired for
        # samples
        if self.rtype == 'sample' and self.noresults:
            no_res = self.session.query(NoResult).filter(
                NoResult.sample_id.in_(self.rids))
            for seq in no_res:
                data = self.get_selected_data(seq)
                yield self.writer.format_sequence(data, pseudo_seq=True)
