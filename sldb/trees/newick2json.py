import argparse
import csv
import json
from ete2 import Tree


def _generate_tab_dict(tab_fh):
    seq_info = {}
    for e in csv.DictReader(tab_fh, delimiter='\t'):
        seq_info[e['seqID']] = e

    return seq_info


def _traverse(node, seq_info):
    node_inst = {
        'name': node.name,
        'children': []
    }

    if node.name not in ['NoName', 'germline']:
        node_inst['data'] = {}
        for e in ['subject', 'subset', 'tissue', 'disease']:
            node_inst['data'][e] = seq_info[node.name][e]

    for child in node.children:
        node_inst['children'].append(_traverse(child, seq_info))

    return node_inst


def _collapse_null_branches(node):
    '''
    Collapse 0 distance branches from an ete2 tree
    (`ete2.tree.TreeNode`).
    Applies recursively.

    Returns
    -------
    None
    '''
    # Remove a node's parent iff the node's branch length is 0 and the
    # parent's name is 'NoName', that way we avoid removing named, and
    # thus possibly informative, nodes.
    if (node.dist == 0 and hasattr(node.up, 'name') and
        node.up.name == 'NoName'):
        parent = node.up
        grandparent = parent.up
        # node.detach()
        for child in [child for child in parent.get_children() if
                child != node]:
            # Remove all children except the node whose branch length is
            # 0.
            # Ensure that the non-zero distances are preserved.
            dist = child.dist
            node.add_child(child.detach(), dist=dist)
        # Make sure that the children are all detached from the parent.
        # This can be removed once tested
        # print(parent.get_children())
        # assert(len(parent.get_children()) == 0)
        # Remove the empty parent, connecting the node to the
        # grandparent with the branch length preserved.
        parent.delete(preserve_branch_length=True)
    # Recurse through all children
    for child in node.get_children():
        _collapse_null_branches(child)
    return None


def convert(tab_fh, newick_fh):
    t = Tree(newick_fh.read())
    t.set_outgroup(t.search_nodes(name='germline')[0])

    seq_info = _generate_tab_dict(tab_fh)
    _collapse_null_branches(t.get_tree_root())
    return _traverse(t.get_tree_root(), seq_info)


def run_newick2json():
    parser = argparse.ArgumentParser(description='Newick to JSON')
    parser.add_argument('tab_file')
    parser.add_argument('newick_file')

    args = parser.parse_args()

    with open(args.tab_file, 'rU') as tab_fh:
        with open(args.newick_file, 'rU') as newick_fh:
            print json.dumps(convert(tab_fh, newick_fh))
