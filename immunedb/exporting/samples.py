from sqlalchemy import func

from immunedb.common.models import CloneStats, Sample, SampleMetadata
from immunedb.exporting.tsv_writer import StreamingTSV, write_tsv
from immunedb.util.log import logger


class Passthrough:
    def __getattr__(self, attr):
        return 0


def get_samples(session, for_update=False):
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
    for sample in session.query(Sample).order_by(Sample.name):
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


def write_samples(session, args):
    logger.info('Exporting samples')
    write_tsv('samples.tsv', get_samples, session, args.for_update)
