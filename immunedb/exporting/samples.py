from sqlalchemy import func

from immunedb.common.models import CloneStats, Sample, SampleMetadata
from immunedb.exporting.tsv_writer import StreamingTSV
from immunedb.exporting.writer import ExportWriter
from immunedb.util.log import logger


class Passthrough:
    def __getattr__(self, attr):
        return 0


def get_samples(session, for_update=False, sample_ids=None):
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

    if for_update:
        fields = ['name', 'new_name']
    else:
        fields = ['id', 'name', 'subject', 'input_sequences', 'identified',
                  'in_frame', 'stops', 'functional', 'clones']
    fields.extend(meta)
    writer = StreamingTSV(fields)
    yield writer.writeheader()
    samples = session.query(Sample)
    if sample_ids:
        samples = samples.filter(Sample.id.in_(sample_ids))
    for sample in samples.order_by(Sample.name):
        row = {
            'id': sample.id,
            'name': sample.name,
            'new_name': sample.name
        }
        stats = sample.stats if sample.stats else Passthrough()
        if not for_update:
            row.update({
                'subject': sample.subject.identifier,
                'input_sequences': stats.sequence_cnt +
                stats.no_result_cnt,
                'identified': stats.sequence_cnt,
                'in_frame': stats.in_frame_cnt,
                'stops': stats.stop_cnt,
                'functional': stats.functional_cnt,
                'clones': clone_cnts.get(sample.id, 0)
            })

        row.update(sample.metadata_dict)
        yield writer.writerow(row)


def write_samples(session, sample_ids=None, for_update=False, zipped=False,
                  **kwargs):
    logger.info('Exporting samples')
    with ExportWriter(zipped) as fh:
        fh.set_filename('samples.tsv')
        fh.write(get_samples(session, for_update, sample_ids))
        return fh.get_zip_value()
