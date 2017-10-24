import csv
import os

from sqlalchemy.sql import exists

from immunedb.common.models import Sample, Sequence
from immunedb.util.log import logger

REQUIRED_FIELDS = ('file_name', 'study_name', 'sample_name', 'date', 'subject')
OPTIONAL_FIELDS = ('subset', 'tissue', 'disease', 'lab', 'experimenter',
                   'ig_class', 'v_primer', 'j_primer')


class MetadataException(Exception):
    pass


def check_populated(row):
    missing = [f for f in REQUIRED_FIELDS if f not in row]
    if len(missing) > 0:
        if 'sample_name' in missing:
            raise MetadataException('Sample name cannot be blank')
        raise MetadataException(
            'Fields {} cannot be blank for sample {}'.format(
                ','.join(missing), row['sample_name']))


def parse_metadata(session, fh, warn_existing, path):
    reader = csv.DictReader(fh, delimiter='\t')
    provided_fields = set(reader.fieldnames)
    missing_fields = set(REQUIRED_FIELDS) - provided_fields
    if len(missing_fields) > 0:
        raise MetadataException(
            'Metadata is missing the following headers: {}'.format(
                ','.join(missing_fields)))
    unknown_fields = provided_fields - (
        set(REQUIRED_FIELDS).union(set(OPTIONAL_FIELDS)))
    if len(unknown_fields) > 0:
        logger.warning('Ignoring unknown headers in metadata: {}'.format(
            ','.join(unknown_fields)))

    metadata = {}
    for row in reader:
        row = {k: v for k, v in row.iteritems()
               if v is not None and len(v) > 0}
        if len(row) == 0:
            continue
        check_populated(row)
        # Check if the sample name is unique
        if row['sample_name'] in metadata:
            logger.error('Duplicate sample name {} in metadata.'.format(
                row['sample_name']))

        # Check if a sample with the same name is in the database
        sample_in_db = session.query(Sample).filter(
            Sample.name == row['sample_name'],
            exists().where(
                Sequence.sample_id == Sample.id
            )).first()
        if sample_in_db:
            message = 'Sample {} already exists. {}'.format(
                row['sample_name'],
                'Skipping.' if warn_existing else 'Cannot continue.'
            )
            if warn_existing:
                logger.warning(message)
                continue
            else:
                raise MetadataException(message)

        # Check if specified file exists
        if not os.path.isfile(os.path.join(path, row['file_name'])):
            raise MetadataException(
                'File {} for sample {} does not exist'.format(
                    row['file_name'], row['sample_name']))

        metadata[row['sample_name']] = row

    return metadata
