from collections import Counter
import io
import zipfile

from immunedb.common.models import Clone, CloneStats, Sample
from immunedb.exporting.tsv_writer import StreamingTSV, write_tsv
from immunedb.util.log import logger
from immunedb.util.lookups import aas_from_nts


def get_clone_features(session):
    return {c.id: (c.v_gene, c.j_gene, c.cdr3_nt)
            for c in session.query(Clone.id, Clone.v_gene, Clone.j_gene,
                                   Clone.cdr3_nt)}


def get_sample_vdjtools(session, sample, min_clone_size, clone_features):
    writer = StreamingTSV(['count', 'freq', 'cdr3nt', 'cdr3aa', 'v', 'd', 'j'])

    sample_clones = Counter()
    stats = session.query(
        CloneStats.clone_id, CloneStats.total_cnt
    ).filter(
        CloneStats.sample_id == sample.id
    )

    for stat in stats:
        sample_clones[clone_features[stat.clone_id]] += stat.total_cnt

    total = sum(sample_clones.values())
    yield writer.writeheader()
    for key in sorted(sample_clones, key=sample_clones.get, reverse=True):
        counts = sample_clones[key]
        if counts < min_clone_size:
            continue
        v, j, cdr3_nt = key
        yield writer.writerow({
            'count': counts,
            'freq': counts / total,
            'cdr3nt': cdr3_nt,
            'cdr3aa': aas_from_nts(cdr3_nt),
            'v': v,
            'd': '.',
            'j': j,
        })


def write_vdjtools(session, args):
    clone_features = get_clone_features(session)
    for sample in session.query(Sample):
        logger.info('Exporting VDJTools format for sample {}'.format(
            sample.name))
        write_tsv(
            '{}.sample.txt'.format(sample.name),
            get_sample_vdjtools,
            session, sample, args.min_clone_size, clone_features
        )


def get_vjdtools_as_zip(session, min_clone_size=0):
    clone_features = get_clone_features(session)

    out = io.BytesIO()
    with zipfile.ZipFile(out, 'w') as fh:
        for sample in session.query(Sample):
            output = ''.join(get_sample_vdjtools(
                session, sample, min_clone_size, clone_features
            ))
            fh.writestr('{}.sample.txt'.format(sample.name), output)
    yield out.getvalue()
