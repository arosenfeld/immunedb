from sqlalchemy import func

from sldb.common.models import Clone, CloneStats, Sequence


def sample_stats(session, clone_id, force):
    for cstat in session.query(
            Sequence.sample_id,
            func.count(Sequence.seq_id).label('unique'),
            func.sum(Sequence.copy_number).label('total'))\
            .filter(Sequence.clone_id == clone_id)\
            .group_by(Sequence.sample_id):

        existing = session.query(CloneStats).filter(
            CloneStats.clone_id == clone_id,
            CloneStats.sample_id == cstat.sample_id).first() is not None
        if existing:
            if force:
                session.delete(CloneStats).filter(
                    CloneStats.clone_id == clone_id)
            else:
                print '\tSkipping clone {} for sample {}.'.format(
                    clone_id, cstat.sample_id)
                continue

        session.add(CloneStats(
            clone_id=clone_id,
            sample_id=cstat.sample_id,
            unique_cnt=cstat.unique,
            total_cnt=cstat.total))


def run_clone_stats(session, args):
    for i, cid in enumerate(session.query(Clone.id)):
        sample_stats(session, cid.id, args.force)

        if i > 0 and i % 1000 == 0:
            print 'Committing {}'.format(i)
            session.commit()
