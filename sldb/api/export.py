import sldb.util.funcs as funcs
from sldb.common.models import DuplicateSequence, NoResult, SequenceMapping

class SequenceExport(object):
    '''A class to handle exporting sequences in various formats.

    :param Session session: A handle to a database session
    :param str eformat: The export format to use.  Currently "tab", "orig", and
        "clip" for tab-delimited, FASTA, FASTA with filled in germlines, and
        FASTA in CLIP format respectively
    :param str rtype: The type of record to filter the query on.  Currently
        either "sample" or "clone"
    :param int-array rids: A list of IDs of ``rtype`` to export
    :param str-array selected_fields: A list of fields to export
    :param bool duplicates: If duplicate sequences should be included in the
        output
    :param bool noresults: If sequences which could not be assigned a V or J
        should be included in the output

    '''
    _export_fields = [
        'seq_id',
        'identity_seq_id',
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

        'alignment',
        'levenshtein_dist',

        'num_gaps',
        'pad_length',
        
        'v_match',
        'v_length',
        'j_match',
        'j_length',

        'pre_cdr3_match',
        'pre_cdr3_length',
        'post_cdr3_match',
        'post_cdr3_length',

        'in_frame',
        'functional',
        'stop',
        'copy_number',
        
        'sequence',
        ('sequence_filled', lambda seq: seq.identity_seq.sequence_replaced),
        ('germline', lambda seq: seq.identity_seq.germline),

        ('v_call', lambda seq: seq.identity_seq.v_call),
        ('j_call', lambda seq: seq.identity_seq.j_call),
        ('cdr3_nt', lambda seq: seq.identity_seq.junction_nt),
        ('cdr3_aa', lambda seq: seq.identity_seq.junction_aa),
        ('gap_method', lambda seq: seq.identity_seq.gap_method),

        'clone_id',
        ('clone_group_id', lambda seq: seq.clone.group.id),
        ('clone_cdr3_nt', lambda seq: seq.clone.cdr3_nt),
        ('clone_cdr3_aa', lambda seq: seq.clone.group.cdr3_aa),
        ('clone_json_tree', lambda seq: seq.clone.tree),
    ]

    def __init__(self, session, eformat, rtype, rids, selected_fields,
                 duplicates, noresults):
        self.session = session
        self.eformat = eformat
        self.rtype = rtype
        self.rids = rids
        self.selected_fields = selected_fields
        self.duplicates = duplicates
        self.noresults = noresults

    @property
    def export_fields(self):
        return self._export_fields

    def _name_and_field(self, e):
        if type(e) == str:
            return e, lambda seq: getattr(seq, e)
        return e

    def _fasta_entry(self, seq_id, keys, sequence):
        return '>{}{}{}\n{}\n'.format(
            seq_id, 
            '|' if len(keys) > 0 else '',
            '|'.join(map(lambda (k, v): '{}={}'.format(k, v), keys)),
            sequence)

    def _tab_entry(self, data):
        return '{}\n'.format('\t'.join(map(lambda s: str(s[1]), data)))

    def _get_selected_data(self, seq, **overrides):
        data = []
        for field in self.export_fields:
            n, f = self._name_and_field(field)
            if n in self.selected_fields:
                try:
                    if n in overrides:
                        data.append((n, overrides[n]))
                    elif 'clone' not in n or (seq.clone is not None):
                        data.append((n, f(seq)))
                    else:
                        data.append((n, 'None'))
                except:
                    data.append((n, 'None'))
        return data

    def get_data(self):
        # Get all the sequences matching the request
        seqs = self.session.query(SequenceMapping).filter(
            getattr(
                SequenceMapping, '{}_id'.format(self.rtype)
            ).in_(self.rids))

        # If it's a CLIP file, order by identity_seq_id to minimize
        # repetition of germline entries
        if self.eformat == 'clip':
            seqs = seqs.order_by(SequenceMapping.clone_id)
        else:
            # This probably isn't necessary but is a safe guard since we
            # page_query is used and order may not be deterministic on all
            # storage engines
            seqs = seqs.order_by(SequenceMapping.seq_id)

            # If it's a tab file, add the headers based on selected fields
            if self.eformat == 'tab':
                headers = []
                for field in self.export_fields:
                    n, f = self._name_and_field(field)
                    if n in self.selected_fields:
                        headers.append(n)
                yield '{}\n'.format('\t'.join(headers))

        # For CLIP files to check if the germline needs to be output
        last_iden = None
        # TODO: Change this to .yield_per?
        for seq in funcs.page_query(seqs):
            # Get the selected data for the sequence
            data = self._get_selected_data(seq)
            if self.eformat == 'fill':
                seq_nts = seq.identity_seq.sequence_replaced
            else:
                seq_nts = seq.sequence

            # If it's a tab file, just output the row
            if self.eformat == 'tab':
                yield self._tab_entry(data)
            else:
                # If it's a CLIP file and there has been a germline change
                # output the new germline with metadata in the header
                if self.eformat == 'clip' and last_iden != seq.identity_seq_id:
                    last_iden = seq.identity_seq_id
                    yield self._fasta_entry(
                        '>Germline',
                        (('v_gene', seq.identity_seq.v_call),
                         ('j_gene', seq.identity_seq.j_call),
                         ('cdr3_aa', seq.identity_seq.junction_aa),
                         ('cdr3_len', seq.identity_seq.junction_num_nts)),
                        seq.identity_seq.germline)

                # Output the FASTA row
                yield self._fasta_entry(seq.seq_id, data, seq_nts)

            # Regardless of output type, output duplicates if desired
            if self.duplicates and seq.copy_number > 1:
                for dup in self.session.query(DuplicateSequence).filter(
                        DuplicateSequence.duplicate_seq == seq):
                    # Duplicates have the same data as their original sequence,
                    # except the seq_id, copy_number is set to 0, and the
                    # duplicate_of_seq_id is set to the original sequence's
                    data = self._get_selected_data(
                        seq,
                        seq_id=dup.seq_id,
                        copy_number=0,
                        duplicate_of_seq_id=seq.seq_id)

                    if self.eformat == 'tab':
                        yield self._tab_entry(data)
                    else:
                        # Output the FASTA row.  Note we don't need to re-output
                        # a germline since these are duplicates and will share
                        # the same as its primary sequence
                        yield self._fasta_entry(dup.seq_id, data, seq_nts)

        # Regardless of the output type, output noresults if desired for samples
        if self.rtype == 'sample' and self.noresults:
            no_res = self.session.query(NoResult).filter(
                NoResult.sample_id.in_(self.rids))
            for seq in no_res:
                data = self._get_selected_data(seq)
                if self.eformat == 'tab':
                    yield self._tab_entry(data)
                else:
                    yield self._fasta_entry(seq.seq_id, data, seq.sequence)
