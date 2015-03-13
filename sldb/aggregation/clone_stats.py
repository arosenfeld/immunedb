import json

from sqlalchemy import func

from sldb.common.models import Clone, CloneStats, Sequence
from sldb.common.mutations import CloneMutations


def clone_stats(session, clone_id, force):
    mutations = {}
    for cstat in session.query(
            Sequence.sample_id,
            func.count(Sequence.seq_id).label('unique'),
            func.sum(Sequence.copy_number).label('total'))\
            .filter(Sequence.clone_id == clone_id)\
            .group_by(Sequence.sample_id):

        existing = session.query(CloneStats).filter(
            CloneStats.clone_id == clone_id,
            CloneStats.sample_id == cstat.sample_id).first()
        if existing is not None:
            if force:
                session.query(CloneStats).filter(
                    CloneStats.clone_id == clone_id,
                    CloneStats.sample_id == sample_id).delete()
            else:
                continue

        if clone_id not in mutations:
            mutations[clone_id] = CloneMutations(
                session,
                session.query(Clone).filter(Clone.id == clone_id).first()
            ).calculate(commit_seqs=True)

        sample_muts = mutations[clone_id][cstat.sample_id]

        session.add(CloneStats(
            clone_id=clone_id,
            sample_id=cstat.sample_id,
            unique_cnt=cstat.unique,
            total_cnt=cstat.total,
            mutations=json.dumps(sample_muts)))


def run_clone_stats(session, args):
    if args.clone_ids is None:
        clones = map(lambda c: c.id, session.query(Clone.id).all())
    else:
        clones = args.clone_ids

    for i, cid in enumerate(clones):
        clone_stats(session, cid, args.force)

        if i > 0 and i % 1000 == 0:
            print 'Committing {}'.format(i)
            session.commit()
    session.commit()
