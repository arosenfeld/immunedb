import shlex
from subprocess import Popen, PIPE

import ete3

import immunedb.common.config as config
from immunedb.common.models import Clone
import immunedb.common.modification_log as mod_log
from immunedb.trees import instantiate_node, LineageWorker
import immunedb.util.concurrent as concurrent
from immunedb.util.log import logger


def get_newick(fasta_input):
    proc = Popen(shlex.split('clearcut --alignment -q --DNA -N -r'),
                 stdin=PIPE, stdout=PIPE, stderr=PIPE, encoding='utf8')
    return proc.communicate(input=fasta_input)[0]


def minimize_tree(tree):
    tree.set_outgroup('germline')
    tree.search_nodes(name='germline')[0].delete()
    first = True
    while True:
        push_common_mutations_up(tree, first)
        remove_parent_mutations(tree)
        rem_null = remove_null_nodes(tree)
        moved = check_supersets(tree)
        if not rem_null and not moved and not are_null_nodes(tree):
            break
        first = False
    return tree


def push_common_mutations_up(tree, first):
    if tree.is_leaf():
        return tree.mutations

    common_muts = None
    for child in tree.children:
        child_muts = push_common_mutations_up(child, first)
        if common_muts is None:
            common_muts = child_muts
        else:
            common_muts = common_muts.intersection(child_muts)

    if len(tree.seq_ids) == 0:
        if first:
            tree.mutations = common_muts or tree.mutations
        else:
            tree.mutations = common_muts.union(tree.mutations)

    return tree.mutations


def remove_parent_mutations(tree):
    for node in tree.traverse(strategy='postorder'):
        if node.up is not None and node.up.name != 'germline':
            node.mutations = node.mutations.difference(node.up.mutations)


def remove_null_nodes(tree):
    removed = False
    for node in tree.traverse():
        if node.up is not None and len(node.mutations) == 0:
            node.up.seq_ids.update(node.seq_ids)
            node.up.copy_number += node.copy_number
            node.delete(prevent_nondicotomic=False)
            removed = True
    return removed


def check_supersets(tree):
    if tree.is_leaf():
        return False

    moved = False
    for c1 in tree.children:
        for c2 in tree.children:
            if c1 == c2:
                continue
            if c1.mutations.issubset(c2.mutations):
                c1.detach()
                c2.add_child(c1)
                moved = True
            elif c2.mutations.issubset(c1.mutations):
                c2.detach()
                c1.add_child(c2)
                moved = True

            overlap = c1.mutations.intersection(c2.mutations)
            if len(overlap) > 0:
                c1.detach()
                c2.detach()
                intermediate = instantiate_node(ete3.Tree(name='NoName'))
                intermediate.mutations = overlap
                intermediate.add_child(c1)
                intermediate.add_child(c2)
                tree.add_child(intermediate)
            moved = moved or check_supersets(c1)

    return moved


def are_null_nodes(tree):
    for node in tree.traverse():
        if node.up is not None and len(node.mutations) == 0:
            return True
    return False


def run_clearcut(session, args):
    if args.clone_ids is not None:
        clones = session.query(Clone.id).filter(
            Clone.id.in_(args.clone_ids))
    else:
        if args.subject_ids is not None:
            clones = session.query(Clone.id).filter(
                Clone.subject_id.in_(args.subject_ids))
        else:
            clones = session.query(Clone.id)

    if not args.force:
        clones = clones.filter(Clone.tree.is_(None))
    clones = [c.id for c in clones]
    mod_log.make_mod('clone_tree', session=session, commit=True,
                     info=vars(args))

    tasks = concurrent.TaskQueue()

    logger.info('Creating task queue for clones')
    for clone_id in clones:
        tasks.add_task(clone_id)

    for _ in range(0, args.nproc):
        session = config.init_db(args.db_config)
        tasks.add_worker(LineageWorker(
            session, get_newick,
            args.min_mut_copies, args.min_mut_samples,
            args.min_seq_copies,
            args.min_seq_samples,
            args.exclude_stops,
            post_tree_hook=minimize_tree))

    tasks.start()
