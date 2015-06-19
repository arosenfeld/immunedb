from sldb.common.models import Clone, CloneStats, Sample
from sldb.exporting.export import Exporter
from sldb.util.nested_writer import NestedCSVWriter


class CloneExport(Exporter):
    _forced_fields = ['clone_id', 'sample_id', 'unique_sequences',
                      'total_sequences']

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



