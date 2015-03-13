import json

from sqlalchemy import distinct, func

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

        # First encountering a new clone
        if clone_id not in mutations:
            # Get the per-sample and total mutations.  Commit the mutations for
            # the sequences.
            mutations[clone_id], all_muts = CloneMutations(
                session,
                session.query(Clone).filter(Clone.id == clone_id).first()
            ).calculate(commit_seqs=True)
            # Get the total and unique counts for the entire clone
            counts = session.query(
                func.count(
                    distinct(Sequence.sequence_replaced)
                ).label('unique'),
                func.sum(Sequence.copy_number).label('total')
            ).filter(Sequence.clone_id == clone_id).first()

            # Add the statistics for the whole clone, denoted with a 0 in the
            # sample_id field
            session.add(CloneStats(
                clone_id=clone_id,
                sample_id=0,
                unique_cnt=counts.unique,
                total_cnt=counts.total,
                mutations=json.dumps({
                    'regions': all_muts.region_muts,
                    'positions': all_muts.position_muts
                })
            ))

        sample_muts = mutations[clone_id][cstat.sample_id]

        session.add(CloneStats(
            clone_id=clone_id,
            sample_id=cstat.sample_id,
            unique_cnt=cstat.unique,
            total_cnt=cstat.total,
            mutations=json.dumps({
                'regions': sample_muts.region_muts,
                'positions': sample_muts.position_muts
            })
        ))


def run_clone_stats(session, args):
    if args.clone_ids is not None:
        clones = args.clone_ids
    elif args.subjects is not None:
        clones = map(lambda c: c.id, session.query(Clone.id).filter(
            Clone.subject_id.in_(args.subjects)).all())
    else:
        clones = map(lambda c: c.id, session.query(Clone.id).all())

    for i, cid in enumerate(clones):
        clone_stats(session, cid, args.force)

        if i > 0 and i % 1000 == 0:
            print 'Committing {}'.format(i)
            session.commit()
    session.commit()
