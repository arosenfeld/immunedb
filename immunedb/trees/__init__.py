import json
import shlex

from subprocess import Popen, PIPE

import ete2

from immunedb.common.models import Sequence, SequenceCollapse


class PhylogeneticTree(object):
    def __init__(self, germline_sequence, sequences, min_mut_occurrence=1,
                 min_mut_samples=1):
        self.germline_sequence = germline_sequence
        self.sequences = sequences
        self.min_mut_occurrence = min_mut_occurrence
        self.min_mut_samples = min_mut_samples

    def run(self, session, clearcut_path):
        fasta, removed_muts = get_fasta_input(
            self.germline_sequence, self.sequences,
            min_mut_occurrence=self.min_mut_occurrence,
            min_mut_samples=self.min_mut_samples)

        newick = get_newick(fasta, clearcut_path)
        tree = populate_tree(session, newick, self.germline_sequence,
                             removed_muts)
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

        self.tree = tree
        return self.tree


def populate_tree(session, newick, germline_seq, removed_muts):
    tree = ete2.Tree(newick)
    for node in tree.traverse():
        if node.name not in ('NoName', 'germline', ''):
            seq = session.query(Sequence).filter(
                Sequence.ai == node.name
            ).first()
            seq_ids = {}
            for collapsed_seq in get_seqs_collapsed_to(session, seq):
                seq_ids[collapsed_seq.seq_id] = {
                    'ai': collapsed_seq.ai,
                    'tissue': collapsed_seq.sample.tissue,
                    'subset': collapsed_seq.sample.subset,
                    'ig_class': collapsed_seq.sample.ig_class,
                    'copy_number': collapsed_seq.copy_number,
                    'sample_name': collapsed_seq.sample.name,
                    'sample_id': collapsed_seq.sample.id
                }

            node.name = seq.seq_id
            node.add_feature('seq_ids', seq_ids)
            node.add_feature('copy_number', sum(
                [s['copy_number'] for s in seq_ids.values()]
            ))
            modified_seq = remove_muts(seq.sequence,
                                       removed_muts, germline_seq)
            node.add_feature('mutations', get_mutations(
                germline_seq, modified_seq,
                map(int, json.loads(seq.mutations_from_clone).keys())))
        else:
            node = instantiate_node(node)

    return tree


def get_fasta_input(germline_seq, sequences, min_mut_occurrence=1,
                    min_mut_samples=1):
    seqs = {}
    mut_counts = {}
    for seq in sequences:
        seqs[seq.ai] = seq.clone_sequence
        if seq.mutations_from_clone is None:
            raise Exception(
                'Mutation information not available for sequence '
                '{}. Skipping this clone tree. Was '
                'immunedb_clone_stats run?'.format(seq.seq_id)
            )

        mutations = get_mutations(
            germline_seq, seq.clone_sequence, map(
                int, json.loads(seq.mutations_from_clone).keys())
        )
        for mut in mutations:
            if mut not in mut_counts:
                mut_counts[mut] = {'count': 0, 'samples': set([])}
            mut_counts[mut]['count'] += 1
            mut_counts[mut]['samples'].add(seq.sample_id)

    removed_muts = set([])
    for mut, cnts in mut_counts.iteritems():
        if (cnts['count'] < min_mut_occurrence or
                len(cnts['samples']) < min_mut_samples):
            removed_muts.add(mut)

    for seq_id, seq in seqs.iteritems():
        seqs[seq_id] = remove_muts(seq, removed_muts, germline_seq)

    in_data = '>germline\n{}\n'.format(germline_seq)
    for seq_id, seq in seqs.iteritems():
        in_data += '>{}\n{}\n'.format(seq_id, seq)
    return in_data, removed_muts


def get_seqs_collapsed_to(session, seq):
    sample_level_seqs = session.query(SequenceCollapse.seq_ai).filter(
        SequenceCollapse.collapse_to_subject_seq_ai == seq.ai
    )
    for sample_seq in sample_level_seqs:
        yield session.query(
            Sequence
        ).filter(Sequence.ai == sample_seq.seq_ai).one()


def remove_muts(seq, removes, germline_seq):
    for mut in removes:
        loc, _, to = mut
        loc -= 1
        if seq[loc] == to:
            seq = seq[:loc] + germline_seq[loc] + seq[loc + 1:]
    return seq


def get_newick(fasta_input, tree_prog):
    proc = Popen(shlex.split('{} --alignment -q --DNA -N -r'.format(
        tree_prog)), stdin=PIPE, stdout=PIPE, stderr=PIPE)
    return proc.communicate(input=fasta_input)[0]


def get_mutations(s1, s2, positions):
    return set([(i + 1, s1[i], s2[i]) for i in positions if s1[i] != s2[i]])


def instantiate_node(node):
    node.add_feature('seq_ids', {})
    node.add_feature('copy_number', 0)
    node.add_feature('sequence', None)
    node.add_feature('mutations', set([]))

    return node


def get_nested(seqs, key):
    ret = set([])
    for s in seqs.values():
        ret.add(s.get(key, set([])))
    return sorted(list(ret))


def tree_as_dict(tree, root=True):
    node = {
        'data': {
            'seq_ids': tree.seq_ids,
            'copy_number': tree.copy_number,
            'tissues': get_nested(tree.seq_ids, 'tissue'),
            'subsets': get_nested(tree.seq_ids, 'subset'),
            'ig_classes': get_nested(tree.seq_ids, 'ig_class'),
            'mutations': [{
                'pos': mut[0],
                'from': mut[1],
                'to': mut[2],
            } for mut in tree.mutations]
        },
        'children': []
    }
    for child in tree.children:
        node['children'].append(tree_as_dict(child, root=False))

    if not root or (len(tree.mutations) == 0 and len(tree.seq_ids) == 0):
        return node

    return {
        'data': {
            'seq_ids': {},
            'copy_number': 0,
            'tissues': '',
            'subsets': '',
            'ig_classes': '',
            'mutations': [],
        },
        'children': [node]
    }


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
                intermediate = instantiate_node(ete2.Tree(name='NoName'))
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


def get_total_muts(tree):
    muts = 0
    for node in tree.traverse():
        muts += len(node.mutations)
    return muts


def cut_tree(tree, max_muts, d=0):
    max_muts -= len(tree.mutations)
    if tree.is_leaf() or max_muts <= 0:
        return [tree]
    subtrees = []
    for node in tree.children:
        subtrees.extend(cut_tree(node, max_muts, d+1))
    return subtrees


def get_seq_pks(tree):
    seq_pks = set([
        (s['sample_id'], s['ai']) for s in tree.seq_ids.values()
    ])
    for child in tree.children:
        seq_pks.update(get_seq_pks(child))
    return seq_pks
