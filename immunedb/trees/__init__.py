import json

import ete3

from immunedb.common.models import Clone, Sequence, SequenceCollapse
import immunedb.util.concurrent as concurrent
from immunedb.util.log import logger


class LineageWorker(concurrent.Worker):
    def __init__(self, session, newick_generator, min_mut_copies,
                 min_mut_samples, min_seq_copies, min_seq_samples,
                 exclude_stops, post_tree_hook=None):
        self.session = session
        self.newick_generator = newick_generator
        self.min_mut_copies = min_mut_copies
        self.min_mut_samples = min_mut_samples
        self.min_seq_copies = min_seq_copies
        self.min_seq_samples = min_seq_samples
        self.exclude_stops = exclude_stops
        self.post_tree_hook = post_tree_hook

    def get_tree(self, germline, sequences):
        fasta, removed_muts = get_fasta_input(
            germline, sequences,
            self.min_mut_copies, self.min_mut_samples
        )
        newick = self.newick_generator(fasta)
        if not newick:
            return None
        tree = add_tree_metadata(self.session, newick, germline, removed_muts)
        if self.post_tree_hook:
            tree = self.post_tree_hook(tree)
        return tree

    def do_task(self, clone_id):
        clone_inst = self.session.query(Clone).filter(
            Clone.id == clone_id).first()
        if not clone_inst:
            return

        self.info('Running clone {}'.format(clone_inst.id))

        sequences = self.session.query(
            Sequence
        ).join(SequenceCollapse).filter(
            Sequence.clone_id == clone_id,
            SequenceCollapse.copy_number_in_subject >= self.min_seq_copies,
            SequenceCollapse.samples_in_subject >= self.min_seq_samples,
        )

        if self.exclude_stops:
            sequences = sequences.filter(Sequence.stop == 0)

        sequences = sequences.order_by(Sequence.v_length)

        try:
            tree = self.get_tree(clone_inst.consensus_germline, sequences)
            if not tree:
                logger.warning('No sequences to make tree for clone {}'.format(
                    clone_id))
                return
        except Exception as e:
            logger.error('Error running clone {}: {}'.format(clone_id, e))
            return

        for node_id, node in enumerate(tree.traverse()):
            node.add_feature('node_id', node_id)
        final = {
            'info': {
                'min_mut_copies': self.min_mut_copies,
                'min_mut_samples': self.min_mut_samples,
                'min_seq_copies': self.min_seq_copies,
                'min_seq_samples': self.min_seq_samples,
                'exclude_stops': self.exclude_stops
            },
            'tree': tree_as_dict(tree)
        }
        clone_inst.tree = json.dumps(final)
        self.session.add(clone_inst)
        self.session.commit()


def get_fasta_input(germline_seq, sequences, min_mut_copies, min_mut_samples):
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
            mut_counts[mut]['count'] += seq.copy_number
            mut_counts[mut]['samples'].add(seq.sample_id)

    removed_muts = set([])
    for mut, cnts in mut_counts.items():
        if (cnts['count'] < min_mut_copies or
                len(cnts['samples']) < min_mut_samples):
            removed_muts.add(mut)

    for seq_id, seq in seqs.items():
        seqs[seq_id] = remove_muts(seq, removed_muts, germline_seq)

    in_data = '>germline\n{}\n'.format(germline_seq)
    for seq_id, seq in seqs.items():
        in_data += '>{}\n{}\n'.format(seq_id, seq)
    return in_data, removed_muts


def add_tree_metadata(session, newick, germline_seq, removed_muts):
    tree = ete3.Tree(newick)
    for node in tree.traverse():
        if node.name not in ('NoName', 'germline', ''):
            seq = session.query(Sequence).filter(
                Sequence.ai == node.name
            ).first()
            seq_ids = {}
            for collapsed_seq in get_seqs_collapsed_to(session, seq):
                seq_ids[collapsed_seq.seq_id] = {
                    'ai': collapsed_seq.ai,
                    'copy_number': collapsed_seq.copy_number,
                    'sample_name': collapsed_seq.sample.name,
                    'sample_id': collapsed_seq.sample.id,
                    'metadata': collapsed_seq.sample.metadata_dict
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


def get_mutations(s1, s2, positions):
    return set([(i + 1, s1[i], s2[i]) for i in positions if s1[i] != s2[i]])


def instantiate_node(node):
    node.add_feature('seq_ids', {})
    node.add_feature('copy_number', 0)
    node.add_feature('sequence', None)
    node.add_feature('mutations', set([]))

    return node


def get_nested(seqs, key):
    return sorted(set(
        [s['metadata'][key] for s in seqs.values() if key in s['metadata']]
    ))


def tree_as_dict(tree, root=True):
    all_meta = set(
        [k for seq in tree.seq_ids.values() for k in seq['metadata'].keys()]
    )
    node = {
        'data': {
            'seq_ids': tree.seq_ids,
            'node_id': tree.node_id,
            'copy_number': tree.copy_number,
            'metadata': {
                k: get_nested(tree.seq_ids, k) for k in all_meta
            },
            'mutations': [{
                'pos': mut[0],
                'from': mut[1],
                'to': mut[2],
            } for mut in tree.mutations]
        },
        'children': []
    }
    for i, child in enumerate(tree.children):
        node['children'].append(tree_as_dict(child, root=False))

    if not root or (len(tree.mutations) == 0 and len(tree.seq_ids) == 0):
        return node

    return {
        'data': {
            'seq_ids': {},
            'copy_number': 0,
            'metadata': {},
            'mutations': [],
        },
        'children': [node]
    }


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
