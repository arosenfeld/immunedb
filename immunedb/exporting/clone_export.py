from collections import OrderedDict

from immunedb.common.models import CloneStats
from immunedb.exporting.export import Exporter
from immunedb.util.nested_writer import NestedCSVWriter


def _if_sample(f, val=''):
    def _wrap(stat):
        if stat.sample_id is not None:
            return f(stat)
        return val
    return _wrap


class CloneExport(Exporter):
    allowed_fields = OrderedDict({
        'clone_id': lambda s: s.clone_id,
        'sample_id': _if_sample(lambda s: s.sample_id, 'Total'),
        'unique_sequences': lambda s: s.unique_cnt,
        'total_sequences': lambda s: s.total_cnt,

        'v_gene': lambda s: s.clone.v_gene,
        'j_gene': lambda s: s.clone.j_gene,
        'cdr3_nt': lambda s: s.clone.cdr3_nt,
        'cdr3_aa': lambda s: s.clone.cdr3_aa,
        'cdr3_num_nts': lambda s: s.clone.cdr3_num_nts,
        'functional': lambda s: s.clone.functional,

        'sample_name': _if_sample(lambda s: s.sample.name, 'Total'),
        'subject_id': _if_sample(lambda s: s.sample.subject.id),
        'subject_identifier':
            _if_sample(lambda s: s.sample.subject.identifier),
        'tissue': _if_sample(lambda s: s.sample.tissue),
        'subset': _if_sample(lambda s: s.sample.subset),
        'disease': _if_sample(lambda s: s.sample.disease),
        'ig_class': _if_sample(lambda s: s.sample.ig_class),
        'v_primer': _if_sample(lambda s: s.sample.v_primer),
        'j_primer': _if_sample(lambda s: s.sample.j_primer),
        'lab': _if_sample(lambda s: s.sample.lab),
        'experimenter': _if_sample(lambda s: s.sample.experimenter),
        'date': _if_sample(lambda s: s.sample.date),
        'tree': lambda s: s.clone.tree,
        'parent_id': lambda s: s.clone.parent_id
    })

    def __init__(self, session, rtype, rids, selected_fields):
        super(CloneExport, self).__init__(
            session, rtype, rids, CloneExport.allowed_fields,
            selected_fields)

    def get_data(self):
        csv = NestedCSVWriter(self.selected_fields, streaming=True)
        if self.rtype == 'sample':
            stats = self.session.query(
                CloneStats.clone_id
            ).filter(
                CloneStats.sample_id.in_(self.rids)
            )
            clone_ids = [e.clone_id for e in stats]
        else:
            clone_ids = self.rids

        for clone_id in sorted(set(clone_ids)):
            stats = self.session.query(CloneStats).filter(
                CloneStats.clone_id == clone_id
            ).order_by(CloneStats.sample_id)
            for stat in stats:
                yield csv.add_row(self.get_selected_data(stat))
