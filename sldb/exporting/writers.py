from sldb.common.models import Sequence

class SequenceWriter(object):
    def preprocess_query(self, sequence_query):
        return sequence_query

    def preprocess_headers(self, headers):
        return headers

    def format_sequence(self, data, pseudo_seq=False):
        raise NotImplementedError()


class FASTAWriter(SequenceWriter):
    def __init__(self, replaced_sequences):
        self._replaced_sequences = replaced_sequences

    def format_sequence(self, data):
        return '>{}\n{}\n'.format(
            _get_header(data)
            _get_sequence(data['sequence'], data['germline'])
        )

    def _get_sequence(self, sequence, germline):
        if not self._replaced_sequences:
            return sequence
        return ''.join(
            [g if s == 'N' else s for s, g in zip(sequence, germline)]
        )

    def _get_header(self, metadata):
        header_fields = {
            k: v for k, v in metadata.iteritems() if k not in
                ('seq_id', 'sequence', 'germline', 'quality')
        }
        return ''.join(
            metadata['seq_id'],
            '|' if len(header_fields) > 0 else '',
            '|'.join(map(lambda (k, v): '{}={}'.format(k, v),
                header_fields, iteritems())
            )
        )

class FASTQWriter(FASTAWriter):
    def __init__(self, replaced_sequences):
        super(FASTQWriter, self).__init__(replaced_sequences)

    def format_sequence(self, data, pseudo_seq=False):
        if 'quality' not in data:
            return ''

        return '@{}\n{}\n+{}\n'.format(
            _get_header(data)
            _get_sequence(data['sequence'], data['germline']),
            data['quality']
        )


class CLIPWriter(FASTAWriter):
    def __init__(self, replaced_sequences):
        super(CLIPWriter, self).__init__(replaced_sequences)

    def preprocess_query(self, headers, sequence_query):
        return sequence_query.order_by(Sequence.clone_id)

    def format_sequence(self, data, pseudo_seq=False):
        if data['clone_id'] != self._last_cid and not pseudo_seq:
            self._last_cid = data['clone_id']
            clone_row = super(CLIPWriter, self).format_sequence({
                'seq_id': 'Germline',
                'v_gene': data['v_gene'],
                'j_gene': data['j_gene'],
                'cdr3_aa': data['cdr3_aa'],
                'cdr3_nt': data['cdr3_nt'],
                'cdr3_len': data['cdr3_len'],
                'sequence': data['germline']
                'germline': data['germline']
            })
        else:
            clone_row = ''
        return clone_row + super(CLIPWriter, self).format_sequence(data)

class CSVWriter(SequenceWriter):
    def preprocess_headers(self, headers):
        self._csv = NestedCSVWriter(headers, streaming=True)

    def format_sequence(self, data, pseudo_seq=False):
        return self._csv.add_row(data)
