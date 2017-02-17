from immunedb.common.models import Sequence
from immunedb.util.nested_writer import NestedCSVWriter


class SequenceWriter(object):
    def __init__(self):
        self._required_fields = []

    def set_selected_fields(self, headers):
        self._selected_fields = headers[:]

    def get_required_fields(self, headers):
        for header in self._required_fields:
            if header not in headers:
                headers.append(header)
        return headers

    def preprocess_query(self, sequence_query):
        return sequence_query.order_by(Sequence.ai)

    def format_sequence(self, data, pseudo_seq=False):
        raise NotImplementedError()


class FASTAWriter(SequenceWriter):
    def __init__(self, replaced_sequences):
        self._replaced_sequences = replaced_sequences
        self._required_fields = ['seq_id', 'sequence', 'germline']

    def format_sequence(self, data):
        return '>{}\n{}\n'.format(
            self._get_header(data),
            self._get_sequence(data['sequence'], data['germline'])
        )

    def _get_sequence(self, sequence, germline):
        if not self._replaced_sequences:
            return sequence

        return ''.join(
            [g if s == 'N' else s for s, g in zip(sequence, germline)]
        )

    def _get_header(self, metadata):
        fields = {k: v for k, v in metadata.iteritems() if k in
                  self._selected_fields}
        for field in self._required_fields:
            if field in fields:
                del fields[field]
        return ''.join([
            metadata['seq_id'],
            '|' if len(fields) > 0 else '',
            '|'.join(['{}={}'.format(k, fields[k]) for k in
                     sorted(fields.keys())])
        ])


class FASTQWriter(FASTAWriter):
    def __init__(self, replaced_sequences):
        super(FASTQWriter, self).__init__(replaced_sequences)
        self._required_fields += ['quality']

    def format_sequence(self, data, pseudo_seq=False):
        return '@{}\n{}\n+\n{}\n'.format(
            self._get_header(data),
            self._get_sequence(data['sequence'], data['germline']),
            data['quality'] or ''
        )


class CLIPWriter(FASTAWriter):
    def __init__(self, replaced_sequences):
        super(CLIPWriter, self).__init__(replaced_sequences)
        self._required_fields += ['v_gene', 'j_gene', 'cdr3_aa', 'cdr3_nt',
                                  'cdr3_num_nts']
        self._last_cid = None

    def preprocess_query(self, sequence_query):
        return sequence_query.order_by(Sequence.clone_id, Sequence.ai)

    def format_sequence(self, data, pseudo_seq=False):
        if data['clone_id'] != self._last_cid and not pseudo_seq:
            self._last_cid = data['clone_id']
            clone_row = super(CLIPWriter, self).format_sequence({
                'seq_id': 'Germline',
                'v_gene': data['v_gene'],
                'j_gene': data['j_gene'],
                'cdr3_aa': data['cdr3_aa'],
                'cdr3_nt': data['cdr3_nt'],
                'cdr3_len': data['cdr3_num_nts'],
                'sequence': data['germline'],
                'germline': data['germline']
            })
        else:
            clone_row = ''
        return clone_row + super(CLIPWriter, self).format_sequence(data)


class CSVWriter(SequenceWriter):
    def set_selected_fields(self, headers):
        super(CSVWriter, self).set_selected_fields(headers)
        self._csv = NestedCSVWriter(headers, streaming=True)

    def format_sequence(self, data, pseudo_seq=False):
        return self._csv.add_row(data)
