from sqlalchemy import func

from immunedb.common.models import CloneStats, Sample, SampleMetadata
from immunedb.exporting.tsv_writer import StreamingTSV, write_tsv
from immunedb.util.log import logger


def get_samples(session):
    meta = [
        s.key for s in session.query(SampleMetadata.key).group_by(
            SampleMetadata.key).order_by(SampleMetadata.key)
    ]

    clone_cnts = {s.sample_id: s.clones for s in session.query(
        CloneStats.sample_id,
        func.count(CloneStats.clone_id.distinct()).label('clones')
    ).filter(
        ~CloneStats.sample_id.is_(None)
    ).group_by(CloneStats.sample_id)}

    writer = StreamingTSV(['id', 'name', 'subject', 'input_sequences',
                           'identified', 'in_frame', 'stops', 'functional',
                           'clones'] +
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
            'functional': sample.stats.functional_cnt,
            'clones': clone_cnts[sample.id]
        }
        row.update(sample.metadata_dict)
        yield writer.writerow(row)


def write_samples(session, args):
    logger.info('Exporting samples')
    write_tsv('samples.tsv', get_samples, session)
