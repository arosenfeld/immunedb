from immunedb.common.models import Clone, CloneStats
from immunedb.exporting.tsv_writer import StreamingTSV, write_tsv
from immunedb.util.funcs import yield_limit
from immunedb.util.log import logger


def get_clone_summary(session, include_lineages):
    fields = [
        'clone_id', 'subject', 'v_gene', 'j_gene', 'functional', 'insertions',
        'deletions', 'cdr3_nt', 'cdr3_num_nt', 'cdr3_aa',
        'uniques', 'instances', 'copies', 'germline', 'parent_id',
        'avg_mutations_per_copy'
    ]
    if include_lineages:
        fields.append('lineage')
    writer = StreamingTSV(fields)

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
            'avg_mutations_per_copy': round(
                clone.overall_stats.total_mutations(normalize=True), 2)
        })
        if include_lineages:
            row['lineage'] = clone.tree
        yield writer.writerow(row)


def get_clone_overlap(session):
    writer = StreamingTSV(['clone_id', 'sample', 'uniques', 'copies',
                           'avg_mutations_per_copy'])

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
            'copies': stat.total_cnt,
            'avg_mutations_per_copy': round(
                stat.total_mutations(normalize=True), 2)
        })


def write_clone_summary(session, args):
    logger.info('Exporting clone summary {} lineages'.format(
        'INCLUDING' if args.include_lineages else 'EXCLUDING'))
    write_tsv('clones.tsv', get_clone_summary, session, args.include_lineages)


def write_clone_overlap(session, args):
    logger.info('Exporting clone overlap')
    write_tsv('clone_overlap.tsv', get_clone_overlap, session)
