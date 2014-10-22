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

    for child in node.get_descendants():
        node_inst['children'].append(_traverse(child, seq_info))

    return node_inst


def convert(tab_fh, newick_fh):
    t = Tree(newick_fh.read())
    t.set_outgroup(t.search_nodes(name='germline')[0])

    seq_info = _generate_tab_dict(tab_fh)
    return _traverse(t.get_tree_root(), seq_info)


def run_newick2json():
    parser = argparse.ArgumentParser(description='Newick to JSON')
    parser.add_argument('tab_file')
    parser.add_argument('newick_file')

    args = parser.parse_args()

    with open(args.tab_file, 'rU') as tab_fh:
        with open(args.newick_file, 'rU') as newick_fh:
            print json.dumps(convert(tab_fh, newick_fh))
