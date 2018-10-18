from sqlalchemy.orm import joinedload

from immunedb.common.models import SelectionPressure
from immunedb.exporting.tsv_writer import StreamingTSV
from immunedb.exporting.writer import ExportWriter
from immunedb.util.funcs import yield_limit
from immunedb.util.log import logger


def get_selection(session, filter_type=None, sample_ids=None):
    query = session.query(SelectionPressure).options(
        joinedload(SelectionPressure.clone),
        joinedload(SelectionPressure.sample),
    )
    if filter_type == 'overall':
        query = query.filter(SelectionPressure.sample_id.is_(None))
    elif filter_type == 'samples':
        if sample_ids:
            query.filter(SelectionPressure.sample_id.in_(sample_ids))
        else:
            query = query.filter(~SelectionPressure.sample_id.is_(None))

    base_fields = SelectionPressure.__table__.c.keys()
    base_fields.remove('id')
    base_fields.remove('sample_id')

    writer = StreamingTSV(['sample', 'subject'] + base_fields)
    yield writer.writeheader()

    for sel in yield_limit(query, SelectionPressure.id):
        row = {f: getattr(sel, f) for f in base_fields}
        row['sample'] = sel.sample.name if sel.sample else 'All Samples'
        row['subject'] = sel.clone.subject.identifier
        yield writer.writerow(row)


def write_selection(session, sample_ids=None, filter_type='both', zipped=False,
                    **kwargs):
    logger.info('Exporting selection pressure')
    with ExportWriter(zipped=zipped) as fh:
        fh.set_filename('selection_pressure.tsv')
        fh.write(get_selection(session, filter_type, sample_ids))
        return fh.get_zip_value()
