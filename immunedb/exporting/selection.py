import csv

from sqlalchemy.orm import joinedload

from immunedb.common.models import SelectionPressure
from immunedb.util.funcs import yield_limit
from immunedb.util.log import logger


def write_selection(session, args):
    query = session.query(SelectionPressure).options(
        joinedload(SelectionPressure.clone),
        joinedload(SelectionPressure.sample),
    )
    if args.filter == 'overall':
        query = query.filter(SelectionPressure.sample_id.is_(None))
    elif args.filter == 'samples':
        query = query.filter(~SelectionPressure.sample_id.is_(None))

    base_fields = SelectionPressure.__table__.c.keys()
    base_fields.remove('id')
    base_fields.remove('sample_id')

    writer = csv.DictWriter(open('selection_pressure.tsv', 'w+'),
                            delimiter='\t', fieldnames=['sample', 'subject'] +
                            base_fields)
    writer.writeheader()

    logger.info('Exporting selection pressure')
    for sel in yield_limit(query, SelectionPressure.id):
        row = {f: getattr(sel, f) for f in base_fields}
        row['sample'] = sel.sample.name if sel.sample else 'All Samples'
        row['subject'] = sel.clone.subject.identifier
        writer.writerow(row)
