from immunedb.common.models import Clone, CloneStats
from immunedb.exporting.tsv_writer import StreamingTSV, write_tsv
from immunedb.util.funcs import yield_limit
from immunedb.util.log import logger


def get_clone_info(session):
    writer = StreamingTSV([
        'clone_id', 'subject', 'v_gene', 'j_gene', 'functional', 'insertions',
        'deletions', 'cdr3_nt', 'cdr3_num_nt', 'cdr3_aa',
        'uniques', 'instances', 'copies', 'germline', 'parent_id'
    ])

    yield writer.writeheader()
    for clone in yield_limit(session.query(Clone), Clone.id):
        row = {}
        for field in writer.fieldnames:
            try:
                row[field] = getattr(clone, field)
            except AttributeError:
                pass
        row.update({
            'clone_id': clone.id,
            'subject': clone.subject.identifier,
            'functional': 'T' if clone.functional else 'F',
            'insertions': clone._insertions,
            'deletions': clone._deletions,
            'uniques': clone.overall_unique_cnt,
            'instances': clone.overall_instance_cnt,
            'copies': clone.overall_total_cnt,
        })
        yield writer.writerow(row)


def get_clone_overlap(session):
    writer = StreamingTSV(['clone_id', 'sample', 'uniques', 'copies'])

    stats = session.query(
        CloneStats
    ).filter(
        ~CloneStats.sample_id.is_(None)
    )

    yield writer.writeheader()
    for stat in yield_limit(stats, CloneStats.id):
        yield writer.writerow({
            'clone_id': stat.clone_id,
            'sample': stat.sample.name,
            'uniques': stat.unique_cnt,
            'copies': stat.total_cnt
        })


def write_clone_info(session, args):
    logger.info('Exporting clone information')
    write_tsv('clones.tsv', get_clone_info, session)


def write_clone_overlap(session, args):
    logger.info('Exporting clone overlap')
    write_tsv('clone_overlap.tsv', get_clone_overlap, session)
