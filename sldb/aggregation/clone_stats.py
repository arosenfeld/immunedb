import json

from sqlalchemy import distinct, func

from sldb.common.models import Clone, CloneStats, Sequence
from sldb.common.mutations import CloneMutations


def clone_stats(session, clone_id, force):
    if force:
        session.query(CloneStats).filter(
            CloneStats.clone_id == clone_id).delete()

    # Get the total and unique counts for the entire clone
    counts = session.query(
        func.count(
            distinct(Sequence.sequence_replaced)
        ).label('unique'),
        func.sum(Sequence.copy_number).label('total')
    ).filter(Sequence.clone_id == clone_id).first()

    mutations = {}
    for cstat in session.query(
            Sequence.sample_id,
            func.count(distinct(Sequence.sequence_replaced)).label('unique'),
            func.sum(Sequence.copy_number).label('total'))\
            .filter(Sequence.clone_id == clone_id)\
            .group_by(Sequence.sample_id):

        existing = session.query(CloneStats).filter(
            CloneStats.clone_id == clone_id,
            CloneStats.sample_id == cstat.sample_id).first()

        # First encountering a new clone
        if clone_id not in mutations:
            # Get the per-sample and total mutations.  Commit the mutations for
            # the sequences.
            mutations[clone_id], all_muts = CloneMutations(
                session,
                session.query(Clone).filter(Clone.id == clone_id).first()
            ).calculate(commit_seqs=True)

            # Add the statistics for the whole clone, denoted with a 0 in the
            # sample_id field
            session.add(CloneStats(
                clone_id=clone_id,
                sample_id=0,
                unique_cnt=counts.unique,
                total_cnt=counts.total,
                mutations=json.dumps(all_muts.get_all())
            ))

        sample_muts = mutations[clone_id][cstat.sample_id]

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
            mutations=json.dumps(sample_muts.get_all())
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
