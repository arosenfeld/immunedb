from sqlalchemy.orm import joinedload

from immunedb.common.models import SelectionPressure
from immunedb.exporting.tsv_writer import StreamingTSV, write_tsv
from immunedb.util.funcs import yield_limit
from immunedb.util.log import logger


def get_selection(session, filter_type=None):
    query = session.query(SelectionPressure).options(
        joinedload(SelectionPressure.clone),
        joinedload(SelectionPressure.sample),
    )
    if filter_type == 'overall':
        query = query.filter(SelectionPressure.sample_id.is_(None))
    elif filter_type == 'samples':
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


def write_selection(session, args):
    logger.info('Exporting selection pressure')
    write_tsv('selection_pressure.tsv', get_selection, session, args.filter)
