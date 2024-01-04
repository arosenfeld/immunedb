import json

from sqlalchemy import update

from immunedb.common.models import Clone, Sequence


def get_all_seq_ids(tree):
    ids = set()
    if tree['data']['seq_ids']:
        ids = ids.union(set(tree['data']['seq_ids'].keys()))

    for child in tree['children']:
        ids = ids.union(get_all_seq_ids(child))
    return ids


def find_removals(tree, cutoff, found=None):
    if len(tree['data']['mutations']) >= cutoff:
        return get_all_seq_ids(tree)

    found = set()
    for child in tree['children']:
        found = found.union(find_removals(child, cutoff, found))

    return found


def trim_clones(session, args):
    removals = set()
    clones = (
        session
        .query(Clone.tree)
        .filter(
            Clone.id == 7972
        )
    )
    for clone in clones:
        removals = removals.union(
            find_removals(json.loads(clone.tree)['tree'], args.cutoff)
        )

    session.execute(
        update(Sequence)
        .where(Sequence.seq_id.in_(removals))
        .values(clone_id=None)
    )
    session.commit()
