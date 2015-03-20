import copy
import base64
import json
import shlex
from subprocess import Popen, PIPE, STDOUT

from sqlalchemy import desc, distinct
from sqlalchemy.sql import func

from sldb.common.models import Clone, CloneGroup, Sample, Sequence
import sldb.common.modification_log as mod_log
from sldb.identification.v_genes import VGene

import ete2


def _remove_muts(seq, removes, germline_seq):
    for mut in removes:
        loc, _, to = mut
        loc -= 1
        if seq[loc] == to:
            seq = seq[:loc] + germline_seq[loc] + seq[loc + 1:]
    return seq


def _get_fasta_input(session, germline_seq, clone_id, min_count):
    seqs = {}
    mut_counts = {}
    for i, seq in enumerate(session.query(
            Sequence.seq_id,
            Sequence.sequence_replaced).filter(
                Sequence.clone_id == clone_id)
            .order_by(desc('copy_number'))
            .group_by(Sequence.sequence_replaced)):
        seqs[base64.b64encode(seq.seq_id)] = seq.sequence_replaced

        for mut in _get_mutations(germline_seq, seq.sequence_replaced):
            if mut not in mut_counts:
                mut_counts[mut] = 0
            mut_counts[mut] += 1

    remove_muts = set([])
    for mut, cnt in mut_counts.iteritems():
        if cnt < min_count:
            remove_muts.add(mut)

    for seq_id, seq in seqs.iteritems():
        seqs[seq_id] = _remove_muts(seq, remove_muts, germline_seq)

    in_data = '>germline\n{}\n'.format(germline_seq)
    for seq_id, seq in seqs.iteritems():
        in_data += '>{}\n{}\n'.format(seq_id, seq)
    return in_data, remove_muts


def _get_newick(session, fasta_input, tree_prog):
    proc = Popen(shlex.split('{} --alignment -q --DNA -N -r'.format(
        tree_prog)), stdin=PIPE, stdout=PIPE)
    return proc.communicate(input=fasta_input)[0]


def _get_mutations(s1, s2):
    muts = set([])
    for i, (c1, c2) in enumerate(zip(s1, s2)):
        if (c1 != c2 and c1 != 'N' and c1 != '-' and c2 != '-'):
            muts.add((i + 1, c1, c2))
    return muts

def _instantiate_node(node):
    node.add_feature('seq_ids', [])
    node.add_feature('copy_number', 0)
    node.add_feature('sequence', None)
    node.add_feature('tissues', [])
    node.add_feature('subsets', [])
    node.add_feature('mutations', set([]))

    return node

def _get_tree(session, newick, germline_seq, remove_muts):
    tree = ete2.Tree(newick)
    mut_counts = {}
    for node in tree.traverse():
        if node.name not in ('NoName', 'germline'):
            name = base64.b64decode(node.name)
            seq = session.query(Sequence).filter(
                Sequence.seq_id == name
            ).first()
            node.name = name
            sample_info = session.query(Sequence).filter(
                Sequence.sequence_replaced == seq.sequence_replaced,
                Sequence.sample.has(subject_id=seq.sample.subject_id)
            ).all()
            seq_ids = [seq.seq_id]
            tissues = set(map(lambda s: s.sample.tissue, sample_info))
            subsets = set(map(lambda s: s.sample.subset, sample_info))
            node.add_feature('seq_ids', seq_ids)
            node.add_feature('copy_number',
                    sum(map(lambda seq: seq.copy_number, sample_info)))
            node.add_feature('tissues', map(str, tissues))
            node.add_feature('subsets', map(str, subsets))
            modified_seq = _remove_muts(seq.sequence_replaced, remove_muts,
                                        germline_seq)
            node.add_feature('mutations', _get_mutations(
                             germline_seq, modified_seq))
        else:
            node = _instantiate_node(node)

    return tree


def _get_json(tree, root=True):
    muts = []
    for mut in tree.mutations:
        muts.append({
            'pos': mut[0],
            'from': mut[1],
            'to': mut[2],
        })
    node = {
        'data': {
            'seq_ids': sorted(set(tree.seq_ids)),
            'copy_number': tree.copy_number,
            'tissues': ','.join(sorted(set(tree.tissues))),
            'subsets': ','.join(sorted(set(tree.subsets))),
            'mutations': muts,
        },
        'children': []
    }
    for child in tree.children:
        node['children'].append(_get_json(child, root=False))

    if not root or (len(tree.mutations) == 0 and tree.name != 'NoName'):
        return node

    return {
        'data': {
            'seq_ids': [],
            'copy_number': 0,
            'tissues': '',
            'subsets': '',
            'mutations': [],
        },
        'children': [node]
    }


def _push_common_mutations_up(tree, first):
    if tree.is_leaf():
        return tree.mutations

    common_muts = None
    for child in tree.children:
        child_muts = _push_common_mutations_up(child, first)
        if common_muts is None:
            common_muts = child_muts
        else:
            common_muts = common_muts.intersection(child_muts)

    if first:
        tree.mutations = common_muts or tree.mutations
    else:
        tree.mutations = common_muts.union(tree.mutations)

    return tree.mutations


def _remove_parent_mutations(tree):
    for node in tree.traverse(strategy='postorder'):
        if node.up is not None and node.up.name != 'germline':
            node.mutations = node.mutations.difference(node.up.mutations)


def _remove_null_nodes(tree):
    for node in tree.traverse():
        if node.up is not None and len(node.mutations) == 0:
            node.up.tissues.extend(node.tissues)
            node.up.subsets.extend(node.subsets)
            node.up.seq_ids.extend(node.seq_ids)
            node.up.copy_number += node.copy_number
            node.delete(prevent_nondicotomic=False)


def _check_supersets(tree):
    if tree.is_leaf():
        return False

    moved = False
    dbg = tree.up is None
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
                intermediate = _instantiate_node(ete2.Tree(name='NoName'))
                intermediate.mutations = overlap
                intermediate.add_child(c1)
                intermediate.add_child(c2)
                tree.add_child(intermediate)
            moved = moved or _check_supersets(c1)

    return moved


def _are_null_nodes(tree):
    for node in tree.traverse():
        if node.up is not None and len(node.mutations) == 0:
            return True
    return False


def _get_total_muts(tree):
    muts = 0
    for node in tree.traverse():
        muts += len(node.mutations)
    return muts

def run_nj(session, args):
    if args.clone_ids is None:
        clones = session.query(Clone.id)
        if not args.force:
            clones = clones.filter(Clone.tree.is_(None))
        clones = map(lambda c: c.id, clones)
    else:
        clones = args.clone_ids
    mod_log.make_mod('clone_tree', session=session, commit=True,
                     info=vars(args))

    for clone in clones:
        clone_inst = session.query(Clone).filter(
            Clone.id == clone).first()
        if clone_inst.tree is not None and not args.force:
            print ('Not regenerating tree for clone {}.  Use --force to '
                   'override.').format(clone)
            continue
        germline_seq = session.query(CloneGroup.germline).filter(
            CloneGroup.id == clone_inst.group_id
        ).first().germline

        germline_seq = (germline_seq[:VGene.CDR3_OFFSET] +
                        clone_inst.cdr3_nt +
                        germline_seq[VGene.CDR3_OFFSET +
                                     clone_inst.cdr3_num_nts:])

        fasta, remove_muts = _get_fasta_input(session, germline_seq, clone,
                                              args.min_count)
        newick = _get_newick(session, fasta, args.tree_prog)

        try:
            tree = _get_tree(session, newick, germline_seq, remove_muts)
        except:
            print '[ERROR] Could not get tree for {}'.format(
                clone)
            continue
        tree.set_outgroup('germline')
        tree.search_nodes(name='germline')[0].delete()

        first = True
        while True:
            _push_common_mutations_up(tree, first)
            _remove_parent_mutations(tree)
            _remove_null_nodes(tree)
            moved = _check_supersets(tree)
            if not moved and not _are_null_nodes(tree):
                break
            first = False

        clone_inst.tree = json.dumps(_get_json(tree))
        session.commit()
