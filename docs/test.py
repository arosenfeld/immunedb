from sqlalchemy.sql import func

import immunedb.common.config as config
from immunedb.common.models import (Clone, CloneStats, Sequence,
                                    SequenceCollapse)

session = config.init_db('/home/arosenfeld/configs/develop.json')


for clone in session.query(Clone).filter(Clone.v_gene == 'IGHV3-30'):
    print 'clone {} has AAs {}'.format(clone.id, clone.cdr3_aa)

for stat in session.query(CloneStats).filter(
        CloneStats.clone_id == 31874).order_by(CloneStats.sample_id):
    print 'clone {} has {} unique sequences and {} copies {}'.format(
        stat.clone_id,
        stat.unique_cnt,
        stat.total_cnt,
        ('in sample ' + stat.sample.name) if stat.sample else 'overall')

subject_unique_seqs = session.query(
    func.count(Sequence.seq_id).label('count'),
    Sequence.v_gene
).join(
    SequenceCollapse
).filter(
    Sequence.subject_id == 1,
    ~Sequence.clone_id.is_(None),
    SequenceCollapse.copy_number_in_subject > 0
).group_by(
    Sequence.v_gene
).order_by(
    'count'
)

for seq in subject_unique_seqs:
    print seq.v_gene, seq.count
