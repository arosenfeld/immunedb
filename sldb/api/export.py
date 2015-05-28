import sldb.api.queries as queries
from sldb.common.models import (Clone, CloneStats, DuplicateSequence, NoResult,
                                Sample, Sequence)
from sldb.common.mutations import threshold_mutations
import sldb.util.funcs as funcs
from sldb.util.nested_writer import NestedCSVWriter


class Exporter(object):
    def __init__(self, session, rtype, rids, export_fields, selected_fields):
        self.session = session
        self.rtype = rtype
        self.rids = rids
        self.selected_fields = selected_fields
        self.export_fields = export_fields

    def _name_and_field(self, e):
        """Gets the name and field for export field entry ``e``.

        :param tuple e: Input field to split into name and field.

        :returns: Name and field for export field ``e``
        :rtype: (name, field)

        """
        if type(e) == str:
            return e, lambda r: getattr(r, e)
        return e

    def get_selected_data(self, seq, **overrides):
        """Gets the data specified by ``selected_fields`` for the sequence
        ``seq`` while overriding values in ``overrides`` if they exist.

        :param Sequence seq: The sequence from which to gather fields
        :param kwargs overrides: Fields to override

        :returns: A list of ``(name, value)`` tuples with the selected data
        :rtype: list

        """
        data = {}
        for field in self.export_fields:
            n, f = self._name_and_field(field)
            if n in self.selected_fields:
                try:
                    if n in overrides:
                        data[n] = overrides[n]
                    elif 'clone' not in n or (seq.clone is not None):
                        data[n] = f(seq)
                    else:
                        data[n] = ''
                except:
                    data[n] = ''
        return data


class CloneExport(Exporter):
    _forced_fields = ['clone_id', 'sample_id', 'unique_sequences',
                      'total_sequences']

    _allowed_fields = [
        'clone_id',
        'sample_id',
        ('unique_sequences', lambda s: s.unique_cnt),
        ('total_sequences', lambda s: s.total_cnt),

        ('group_id', lambda s: s.clone.group_id),
        ('v_gene', lambda s: s.clone.v_gene),
        ('j_gene', lambda s: s.clone.j_gene),
        ('cdr3_nt', lambda s: s.clone.cdr3_nt),
        ('cdr3_aa', lambda s: s.clone.group.cdr3_aa),
        ('cdr3_num_nts', lambda s: s.clone.cdr3_num_nts),
        ('functional', lambda s: s.clone.cdr3_num_nts % 3 == 0),

        ('sample_name', lambda s: s.sample.name),
        ('subject_id', lambda s: s.sample.subject.id),
        ('subject_identifier', lambda s: s.sample.subject.identifier),
        ('tissue', lambda s: s.sample.tissue),
        ('subset', lambda s: s.sample.subset),
        ('disease', lambda s: s.sample.disease),
        ('lab', lambda s: s.sample.lab),
        ('experimenter', lambda s: s.sample.experimenter),
        ('date', lambda s: s.sample.date),
    ]

    def __init__(self, session, rtype, rids, selected_fields,
                 include_total_row):
        super(CloneExport, self).__init__(
            session, rtype, rids, CloneExport._allowed_fields,
            selected_fields)
        self._include_total_row = include_total_row

    def get_data(self):
        query = self.session.query(CloneStats).filter(
            getattr(
                CloneStats, '{}_id'.format(self.rtype)
            ).in_(self.rids),
        )
        if len(self.selected_fields) > 0:
            query = query.join(Clone).join(Sample)
        query = query.order_by(CloneStats.clone_id)

        headers = []
        for field in self.export_fields:
            n, f = self._name_and_field(field)
            if n in self.selected_fields:
                headers.append(n)
        csv = NestedCSVWriter(headers, streaming=True)

        last_cid = None
        overall_unique = 0
        overall_total = 0
        for record in query:
            if self._include_total_row:
                if last_cid is None:
                    last_cid = record.clone_id
                elif last_cid != record.clone_id:
                    cnts = self.session.query(
                        CloneStats.unique_cnt,
                        CloneStats.total_cnt
                        ).filter(
                        CloneStats.clone_id == last_cid,
                        CloneStats.sample_id == 0
                        ).first()
                    total_row = {
                        'clone_id': last_cid,
                        'sample_id': 'TOTAL',
                        'unique_sequences': cnts.unique_cnt,
                        'total_sequences': cnts.total_cnt
                    }
                    yield csv.add_raw_row(total_row)
                    last_cid = record.clone_id

            yield csv.add_row(self.get_selected_data(record))


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

        'alignment',
        'probable_indel_or_misalign',

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
        'quality',
        'sequence_filled',
        'germline',

        'v_gene',
        'j_gene',
        ('cdr3_nt', lambda seq: seq.junction_nt),
        ('cdr3_aa', lambda seq: seq.junction_aa),
        ('cdr3_num_nts', lambda seq: seq.junction_num_nts),
        'gap_method',

        'clone_id',
        ('clone_group_id', lambda seq: seq.clone.group.id),
        ('clone_cdr3_nt', lambda seq: seq.clone.cdr3_nt),
        ('clone_cdr3_aa', lambda seq: seq.clone.group.cdr3_aa),
        ('clone_json_tree', lambda seq: seq.clone.tree),

        'copy_number_in_sample',
        'collapse_to_sample_seq_id',

        'copy_number_in_subject',
        'collapse_to_subject_seq_id',
        'collapse_to_subject_sample_id',
    ]

    def __init__(self, session, eformat, rtype, rids, selected_fields,
                 min_copy, duplicates, noresults, level):
        super(SequenceExport, self).__init__(
            session, rtype, rids, SequenceExport._allowed_fields,
            selected_fields)
        self.eformat = eformat
        self.min_copy = min_copy
        self.duplicates = duplicates
        self.noresults = noresults
        self.level = level

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

        :returns: Lines of output for the selected sequences in ``eformat``
        :rtype: str

        """
        # Get all the sequences matching the request
        seqs = self.session.query(Sequence).filter(
            getattr(
                Sequence, '{}_id'.format(self.rtype)
            ).in_(self.rids)
        )
        if self.level == 'all':
            seqs = seqs.filter(
                Sequence.copy_number >= self.min_copy
            )
        else:
            seqs = seqs.filter(
                getattr(Sequence, 'copy_number_in_{}'.format(self.level)) >=
                self.min_copy
            )

        # If it's a CLIP file, order by clone_id to minimize
        # repetition of germline entries
        if self.eformat == 'clip':
            seqs = seqs.order_by(Sequence.clone_id)
        else:
            # This probably isn't necessary but is a safe guard since we
            # page_query is used and order may not be deterministic on all
            # storage engines
            seqs = seqs.order_by(Sequence.seq_id)

            # If it's a csv file, add the headers based on selected fields
            if self.eformat == 'csv':
                headers = []
                for field in self.export_fields:
                    n, f = self._name_and_field(field)
                    if n in self.selected_fields:
                        headers.append(n)
                self._csv = NestedCSVWriter(headers, streaming=True)

        # For CLIP files to check if the germline needs to be output
        last_cid = ''
        # TODO: Change this to .yield_per?
        for seq in funcs.page_query(seqs):
            # Get the selected data for the sequence
            data = self.get_selected_data(seq)
            if self.eformat == 'fill' or self.eformat == 'clip':
                seq_nts = seq.sequence_replaced
            else:
                seq_nts = seq.sequence

            # If it's a csv file, just output the row
            if self.eformat == 'csv':
                yield self._csv.add_row(data)
            else:
                # If it's a CLIP file and there has been a germline change
                # output the new germline with metadata in the header
                if self.eformat == 'clip' and last_cid != seq.clone_id:
                    last_cid = seq.clone_id
                    yield self._fasta_entry(
                        '>Germline', {
                            'v_gene': seq.v_gene,
                            'j_gene': seq.j_gene,
                            'cdr3_aa': seq.junction_aa,
                            'cdr3_nt': seq.junction_nt,
                            'cdr3_len': seq.junction_num_nts
                        }, seq.germline)

                # Output the FASTA row
                yield self._fasta_entry(seq.seq_id, data, seq_nts)

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

                    if self.eformat == 'csv':
                        yield self._csv.add_row(data)
                    else:
                        # Output the FASTA row.  Note we don't need to
                        # re-output a germline since these are duplicates
                        # and will share the same as its primary sequence
                        yield self._fasta_entry(dup.seq_id, data, seq_nts)

        # Regardless of the output type, output noresults if desired for
        # samples
        if self.rtype == 'sample' and self.noresults:
            no_res = self.session.query(NoResult).filter(
                NoResult.sample_id.in_(self.rids))
            for seq in no_res:
                data = self.get_selected_data(seq)
                if self.eformat == 'csv':
                    yield self._csv.add_row(data)
                else:
                    yield self._fasta_entry(seq.seq_id, data, seq.sequence)


class MutationExporter(object):
    def __init__(self, session, clone_ids, limit_sample_ids, thresh_type,
                 thresh_value):
        self._session = session
        self._clone_ids = clone_ids
        self._limit_sample_ids = limit_sample_ids
        self._thresh_type = thresh_type
        self._thresh_value = thresh_value

        headers = ['clone_id', 'sample_id', 'total_seqs', 'region']
        self._numbers = ('unique', 'total')
        self._mtypes = ('synonymous', 'nonsynonymous', 'conservative',
                        'nonconservative', 'unknown')
        for number in self._numbers:
            for mtype in self._mtypes:
                headers.append('{}_{}'.format(mtype, number))

        self._csv = NestedCSVWriter(headers, {
            'nonsynonymous_unique': (
                lambda r: r.get('conservative_unique', 0) +
                r.get('nonconservative_unique', 0)),
            'nonsynonymous_total': (
                lambda r: r.get('conservative_total', 0) +
                r.get('nonconservative_total', 0)),
            'sample_id': lambda r: r['sample_id'] if r['sample_id'] else 'All'
        }, streaming=True)

    def _sample_rows(self, cid, sample_id):
        try:
            all_mutations, total_seqs = queries.get_clone_mutations(
                self._session, cid, sample_id)
        except Exception as ex:
            return

        if self._thresh_type == 'seqs':
            min_seqs = self._thresh_value
        else:
            min_seqs = int(math.ceil(self._thresh_value / 100.0 * total_seqs))

        mutations = threshold_mutations(all_mutations, min_seqs)
        for region, stats in mutations.iteritems():
            row = {
                'clone_id': cid,
                'region': region,
                'sample_id': sample_id,
                'total_seqs': total_seqs,
            }

            for number in self._numbers:
                for mtype, count in stats['counts'][number].iteritems():
                    row['{}_{}'.format(mtype, number)] = count
            self._csv.add_row(row, write_default=True, default=0,
                              write_if_stream=False)

    def get_data(self):
        for cid in self._clone_ids:
            # Write mutations for entire clone
            self._sample_rows(cid, None)
            yield self._csv.get_value()
            # Write per sample mutations
            if self._limit_sample_ids is not None:
                sample_ids = self._limit_sample_ids
            else:
                sample_ids = map(lambda r: r.sample_id, self._session.query(
                    CloneStats.sample_id).filter(
                        CloneStats.clone_id == cid,
                        CloneStats.sample_id != 0).all())
            for sample_id in sample_ids:
                self._sample_rows(cid, sample_id)
                yield self._csv.get_value()
