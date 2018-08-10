from immunedb.common.models import Sample, SampleMetadata
from immunedb.exporting.tsv_writer import StreamingTSV, write_tsv
from immunedb.util.log import logger


def get_samples(session):
    meta = [
        s.key for s in session.query(SampleMetadata.key).group_by(
            SampleMetadata.key).order_by(SampleMetadata.key)
    ]
    writer = StreamingTSV(['id', 'name', 'subject', 'input_sequences',
                           'identified', 'in_frame', 'stops', 'functional'] +
                          meta)
    yield writer.writeheader()
    for sample in session.query(Sample).order_by(Sample.name):
        row = {
            'id': sample.id,
            'name': sample.name,
            'subject': sample.subject.identifier,
            'input_sequences': sample.stats.sequence_cnt,
            'identified': sample.stats.sequence_cnt -
            sample.stats.no_result_cnt,
            'in_frame': sample.stats.in_frame_cnt,
            'stops': sample.stats.stop_cnt,
            'functional': sample.stats.functional_cnt
        }
        row.update(sample.metadata_dict)
        yield writer.writerow(row)


def write_samples(session, args):
    logger.info('Exporting samples')
    write_tsv('samples.tsv', get_samples, session)
