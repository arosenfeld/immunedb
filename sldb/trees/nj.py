import copy
import base64
import json
import shlex
from subprocess import Popen, PIPE, STDOUT

from sqlalchemy import desc, distinct
from sqlalchemy.sql import func

from sldb.common.models import Clone, CloneGroup, Sample, Sequence

import ete2
import networkx as nx


def _get_fasta_input(session, germline_seq, clone_id):
    in_data = '>germline\n{}\n'.format(germline_seq)
    for i, seq in enumerate(session.query(
            Sequence.seq_id,
            Sequence.sequence_replaced).filter(
                Sequence.clone_id == clone_id)
            .order_by(desc('copy_number'))
            .group_by(Sequence.sequence_replaced)):

        in_data += '>{}\n{}\n'.format(base64.b64encode(seq.seq_id),
                                      seq.sequence_replaced)
    return in_data


def _get_newick(session, tree_prog, clone_id):
    clone_inst = session.query(Clone).filter(
        Clone.id == clone_id).first()

    germline_seq = session.query(CloneGroup.germline).filter(
        CloneGroup.id == clone_inst.group_id
    ).first().germline

    stdin = _get_fasta_input(session, germline_seq, clone_id)

    proc = Popen(shlex.split('{} --alignment -q --DNA -N'.format(
        tree_prog)), stdin=PIPE, stdout=PIPE)
    return proc.communicate(input=stdin)[0]


def _get_mutations(s1, s2):
    muts = set([])
    for i, (c1, c2) in enumerate(zip(s1, s2)):
        if (c1 != c2 and c1 != 'N' and c1 != '-'):
            muts.add((i + 1, c1, c2))
    return muts


def _get_tree(session, newick, germline_seq):
    tree = ete2.Tree(newick)
    for node in tree.traverse():
        if node.name not in ('NoName', 'germline'):
            name = base64.b64decode(node.name)
            seq = session.query(Sequence).filter(
                Sequence.seq_id == name
            ).first()
            sample_info = session.query(Sequence).filter(
                Sequence.sequence_replaced == seq.sequence_replaced
            ).all()
            seq_ids = map(lambda s: s.seq_id, sample_info)
            tissues = set(map(lambda s: s.sample.tissue, sample_info))
            subsets = set(map(lambda s: s.sample.subset, sample_info))
            node.add_feature('seq_ids', seq_ids)
            node.add_feature('tissues',
                             ','.join(map(str, tissues)))
            node.add_feature('subsets',
                             ','.join(map(str, subsets)))
            node.add_feature('mutations', _get_mutations(
                germline_seq, seq.sequence_replaced))
        else:
            node.add_feature('seq_ids', [])
            node.add_feature('sequence', None)
            node.add_feature('tissues', [])
            node.add_feature('subsets', [])
            node.add_feature('mutations', set([]))
    return tree


def _get_json(tree):
    muts = []
    for mut in tree.mutations:
        muts.append({
            'pos': mut[0],
            'from': mut[1],
            'to': mut[2],
        })
    node = {
        'data': {
            'seq_ids': tree.seq_ids,
            'copy_number': 5,
            'tissues': tree.tissues,
            'subsets': tree.subsets,
            'mutations': muts,
        },
        'children': []
    }
    for child in tree.children:
        node['children'].append(_get_json(child))

    return node


def _push_common_mutations_up(tree):
    if tree.is_leaf():
        return tree.mutations
    common_muts = None
    rest = False
    for child in tree.children:
        child_muts = copy.copy(_push_common_mutations_up(child))
        if common_muts is None:
            common_muts = child_muts
        else:
            common_muts = common_muts.intersection(child_muts)

    tree.mutations = common_muts
    return common_muts


def _remove_parent_mutations(tree):
    for node in tree.traverse(strategy='postorder'):
        if node.up is not None:
            node.mutations = node.mutations.difference(node.up.mutations)

def _remove_null_nodes(tree):
    for node in tree.traverse():
        if len(node.mutations) == 0:
            node.delete()


def run_nj(session, args):
    if args.clone_ids is None:
        clones = session.query(Clone.id)
        if not args.force:
            clones = clones.filter(Clone.tree.is_(None))
        clones = map(lambda c: c.id, clones)
    else:
        clones = args.clone_ids

    for clone in clones:
        clone_inst = session.query(Clone).filter(
            Clone.id == clone).first()
        if clone_inst.tree is not None and not args.force:
            print ('Not regenerating tree for clone {}.  Use --force to '
                   'override.').format(clone)
            continue
        newick = _get_newick(session, args.tree_prog, clone)

        germline_seq = session.query(CloneGroup.germline).filter(
            CloneGroup.id == clone_inst.group_id
        ).first().germline

        tree = _get_tree(session, newick, germline_seq)
        tree.set_outgroup('germline')
        tree.search_nodes(name='germline')[0].delete()
        _push_common_mutations_up(tree)
        _remove_parent_mutations(tree)
        _remove_null_nodes(tree)

        clone_inst.tree = json.dumps(_get_json(tree))
        session.commit()
