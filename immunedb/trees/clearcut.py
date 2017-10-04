import json
import shlex
from subprocess import Popen, PIPE

import ete2

import immunedb.common.config as config
from immunedb.common.models import Clone, Sequence, SequenceCollapse
import immunedb.common.modification_log as mod_log
import immunedb.util.concurrent as concurrent
from immunedb.util.log import logger


class ClearcutWorker(concurrent.Worker):
    def __init__(self, session, tree_prog, min_count, min_samples,
                 min_seq_copies, exclude_stops):
        self._session = session
        self._tree_prog = tree_prog
        self._min_count = min_count
        self._min_samples = min_samples
        self._min_seq_copies = min_seq_copies
        self._exclude_stops = exclude_stops

    def do_task(self, clone_id):
        clone_inst = self._session.query(Clone).filter(
            Clone.id == clone_id).first()
        if clone_inst is None:
            return

        self.info('Running clone {}'.format(clone_inst.id))

        germline_seq = clone_inst.consensus_germline

        fasta, remove_muts = self._get_fasta_input(germline_seq, clone_inst.id)
        newick = _get_newick(fasta, self._tree_prog)
        try:
            tree = self._get_tree(newick, germline_seq, remove_muts)
        except Exception as e:
            self.error('Error running clone {}: {}'.format(clone_id, e))
            return

        tree.set_outgroup('germline')
        tree.search_nodes(name='germline')[0].delete()

        first = True
        while True:
            _push_common_mutations_up(tree, first)
            _remove_parent_mutations(tree)
            rem_null = _remove_null_nodes(tree)
            moved = _check_supersets(tree)
            if not rem_null and not moved and not _are_null_nodes(tree):
                break
            first = False

        final = {
            'info': {
                'min_count': self._min_count,
                'min_samples': self._min_samples,
                'min_seq_copies': self._min_seq_copies,
                'exclude_stops': self._exclude_stops
            },
            'tree': _get_json(tree)
        }
        clone_inst.tree = json.dumps(final)
        self._session.add(clone_inst)
        self._session.commit()

    def _get_fasta_input(self, germline_seq, clone_id):
        seqs = {}
        mut_counts = {}
        q = self._session.query(
            Sequence
        ).join(SequenceCollapse).filter(
            Sequence.clone_id == clone_id,
            SequenceCollapse.copy_number_in_subject > self._min_seq_copies
        ).order_by(
            Sequence.v_length
        )
        for seq in q:
            if self._exclude_stops and seq.stop:
                continue
            seqs[seq.ai] = seq.clone_sequence
            if seq.mutations_from_clone is None:
                raise Exception(
                    'Mutation information not available for sequence '
                    '{}. Skipping this clone tree. Was '
                    'immunedb_clone_stats run?'.format(seq.seq_id)
                )

            mutations = _get_mutations(
                germline_seq, seq.clone_sequence, map(
                    int, json.loads(seq.mutations_from_clone).keys())
            )
            for mut in mutations:
                if mut not in mut_counts:
                    mut_counts[mut] = {'count': 0, 'samples': set([])}
                mut_counts[mut]['count'] += 1
                mut_counts[mut]['samples'].add(seq.sample_id)

        remove_muts = set([])
        for mut, cnts in mut_counts.iteritems():
            if (cnts['count'] < self._min_count or
                    len(cnts['samples']) < self._min_samples):
                remove_muts.add(mut)

        for seq_id, seq in seqs.iteritems():
            seqs[seq_id] = _remove_muts(seq, remove_muts, germline_seq)

        in_data = '>germline\n{}\n'.format(germline_seq)
        for seq_id, seq in seqs.iteritems():
            in_data += '>{}\n{}\n'.format(seq_id, seq)
        return in_data, remove_muts

    def _get_tree(self, newick, germline_seq, remove_muts):
        tree = ete2.Tree(newick)
        for node in tree.traverse():
            if node.name not in ('NoName', 'germline', ''):
                seq = self._session.query(Sequence).filter(
                    Sequence.ai == node.name
                ).first()
                seq_ids = {}
                for collapsed_seq in _get_seqs_collapsed_to(
                        self._session, seq):
                    seq_ids[collapsed_seq.seq_id] = {
                        'tissue': collapsed_seq.sample.tissue,
                        'subset': collapsed_seq.sample.subset,
                        'ig_class': collapsed_seq.sample.ig_class,
                        'copy_number': collapsed_seq.copy_number,
                        'sample_name': collapsed_seq.sample.name
                    }

                node.name = seq.seq_id
                node.add_feature('seq_ids', seq_ids)
                node.add_feature('copy_number', sum(
                    [s['copy_number'] for s in seq_ids.values()]
                ))
                modified_seq = _remove_muts(seq.sequence, remove_muts,
                                            germline_seq)
                node.add_feature('mutations', _get_mutations(
                    germline_seq, modified_seq,
                    map(int, json.loads(seq.mutations_from_clone).keys())))
            else:
                node = _instantiate_node(node)

        return tree


def _get_seqs_collapsed_to(session, seq):
    sample_level_seqs = session.query(SequenceCollapse.seq_ai).filter(
        SequenceCollapse.collapse_to_subject_seq_ai == seq.ai
    )
    for sample_seq in sample_level_seqs:
        yield session.query(
            Sequence
        ).filter(Sequence.ai == sample_seq.seq_ai).one()


def _remove_muts(seq, removes, germline_seq):
    for mut in removes:
        loc, _, to = mut
        loc -= 1
        if seq[loc] == to:
            seq = seq[:loc] + germline_seq[loc] + seq[loc + 1:]
    return seq


def _get_newick(fasta_input, tree_prog):
    proc = Popen(shlex.split('{} --alignment -q --DNA -N -r'.format(
        tree_prog)), stdin=PIPE, stdout=PIPE, stderr=PIPE)
    return proc.communicate(input=fasta_input)[0]


def _get_mutations(s1, s2, positions):
    return set([(i + 1, s1[i], s2[i]) for i in positions if s1[i] != s2[i]])


def _instantiate_node(node):
    node.add_feature('seq_ids', {})
    node.add_feature('copy_number', 0)
    node.add_feature('sequence', None)
    node.add_feature('mutations', set([]))

    return node


def _get_nested(seqs, key):
    ret = set([])
    for s in seqs.values():
        ret.add(s.get(key, set([])))
    return sorted(list(ret))


def _get_json(tree, root=True):
    node = {
        'data': {
            'seq_ids': tree.seq_ids,
            'copy_number': tree.copy_number,
            'tissues': _get_nested(tree.seq_ids, 'tissue'),
            'subsets': _get_nested(tree.seq_ids, 'subset'),
            'ig_classes': _get_nested(tree.seq_ids, 'ig_class'),
            'mutations': [{
                'pos': mut[0],
                'from': mut[1],
                'to': mut[2],
            } for mut in tree.mutations]
        },
        'children': []
    }
    for child in tree.children:
        node['children'].append(_get_json(child, root=False))

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

    if len(tree.seq_ids) == 0:
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
    removed = False
    for node in tree.traverse():
        if node.up is not None and len(node.mutations) == 0:
            node.up.seq_ids.update(node.seq_ids)
            node.up.copy_number += node.copy_number
            node.delete(prevent_nondicotomic=False)
            removed = True
    return removed


def _check_supersets(tree):
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
        tasks.add_worker(ClearcutWorker(session, args.clearcut_path,
                                        args.min_count, args.min_samples,
                                        args.min_seq_copies,
                                        args.exclude_stops))

    tasks.start()
