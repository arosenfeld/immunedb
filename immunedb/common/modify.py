import csv

from immunedb.common.models import Clone, SampleMetadata

from immunedb.util.log import logger


def update_metadata(session, args):
    with open(args.new_metadata) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        new_meta = {int(l['id']): l for l in reader}
    session.query(SampleMetadata).filter(
        SampleMetadata.sample_id.in_(new_meta)
    ).delete(synchronize_session='fetch')

    base_fields = ['id', 'name', 'subject']
    for sample_id, row in new_meta.items():
        logger.info('Updating sample\'s {} metadata fields'.format(
            row['name']))
        session.add_all([
            SampleMetadata(sample_id=sample_id, key=k, value=v)
            for k, v in row.items() if k not in base_fields and v is not None
        ])

    if session.query(Clone.id).filter(~Clone.tree.is_(None)).count() > 0:
        logger.warning('This database has at least one clonal lineage '
                       'constructed.  All lineages will need to be updated '
                       'to reflect the modified metadata.')
    session.commit()
