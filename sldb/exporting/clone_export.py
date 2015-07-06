from sldb.common.models import Clone, CloneStats, Sample
from sldb.exporting.export import Exporter
from sldb.util.nested_writer import NestedCSVWriter


class CloneExport(Exporter):
    _allowed_fields = [
        'clone_id',
        'sample_id',
        ('unique_sequences', lambda s: s.unique_cnt),
        ('total_sequences', lambda s: s.total_cnt),

        ('v_gene', lambda s: s.clone.v_gene),
        ('j_gene', lambda s: s.clone.j_gene),
        ('cdr3_nt', lambda s: s.clone.cdr3_nt),
        ('cdr3_aa', lambda s: s.clone.cdr3_aa),
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
        ('tree', lambda s: s.clone.tree),
    ]

    # Fields that can be output in total rows
    _total_row_fields = [
        'clone_id', 'sample_id', 'sample_name', 'unique_sequences',
        'total_sequences', 'v_gene', 'j_gene', 'cdr3_nt', 'cdr3_aa',
        'cdr3_num_nts', 'functional'
    ]

    def __init__(self, session, rtype, rids, selected_fields,
                 include_total_row):
        super(CloneExport, self).__init__(
            session, rtype, rids, CloneExport._allowed_fields,
            selected_fields)
        self._include_total_row = include_total_row

    def get_data(self):
        query = self.get_base_query(CloneStats).join(Clone).join(Sample)
        query = query.order_by(CloneStats.clone_id)

        csv = NestedCSVWriter(self.get_headers(), streaming=True)

        last_clone = None
        for record in query:
            if (self._include_total_row and last_clone != record.clone and
                    last_clone is not None):
                yield csv.add_raw_row(self.get_total_row(last_clone))
            last_clone = record.clone
            yield csv.add_row(self.get_selected_data(record))

    def get_total_row(self, clone):
        stats = self.get_selected_data(self.session.query(CloneStats).filter(
            CloneStats.clone == clone,
            CloneStats.sample_id == 0
        ).first())

        total_row = {k: v for k, v in stats.iteritems()
            if k in self._total_row_fields}
        total_row['sample_name'] = 'Total'
        return total_row
