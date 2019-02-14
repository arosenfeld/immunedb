from collections import OrderedDict

from sqlalchemy import desc
from sqlalchemy.sql import text

import dnautils
from immunedb.common.models import (CDR3_OFFSET, Clone, Sequence,
                                    SequenceCollapse, Subject)
import immunedb.common.modification_log as mod_log
import immunedb.common.config as config
from immunedb.trees import cut_tree, get_seq_pks, LineageWorker
import immunedb.trees.clearcut as clearcut
import immunedb.util.concurrent as concurrent
import immunedb.util.funcs as funcs
import immunedb.util.lookups as lookups
from immunedb.util.log import logger


def generate_consensus(session, clone_ids):
    """Generates consensus CDR3s for clones.

    :param Session session: The database session
    :param list clone_ids: The list of clone IDs to assign to groups

    """

    if len(clone_ids) == 0:
        return
    for clone in funcs.periodic_commit(
            session,
            session.query(Clone).filter(Clone.id.in_(clone_ids)),
            interval=1000):
        seqs = session.query(
            Sequence
        ).join(SequenceCollapse).filter(
            Sequence.clone_id == clone.id,
            SequenceCollapse.copy_number_in_subject > 0
        ).all()
        clone.cdr3_nt = funcs.consensus([s.cdr3_nt for s in seqs])
        clone.cdr3_aa = lookups.aas_from_nts(clone.cdr3_nt)

        clone.germline = generate_germline(session, seqs, clone)

    session.commit()


def generate_germline(session, seqs, clone):
    rep_seq = seqs[0]
    cdr3_start_pos = sum(rep_seq.regions[:5])
    germline = rep_seq.germline[:cdr3_start_pos]
    germline += '-' * clone.cdr3_num_nts
    clone.functional = (
        len(germline) % 3 == 0 and
        not lookups.has_stop(germline)
    )

    j_region = rep_seq.germline[cdr3_start_pos + rep_seq.cdr3_num_nts:]
    germline += j_region

    return germline


def push_clone_ids(session):
    session.connection(mapper=Sequence).execute(text('''
        UPDATE
            sequences AS s
        JOIN sequence_collapse AS c
            ON s.sample_id=c.sample_id AND s.ai=c.seq_ai
        JOIN sequences as s2
            ON c.collapse_to_subject_seq_ai=s2.ai
        SET s.clone_id=s2.clone_id
        WHERE s.seq_id!=s2.seq_id
    '''))


def similar_to_all(seq, rest, field, min_similarity):
    """Determines if the string ``seq`` is at least ``min_similarity``
    similar to the list of strings ``rest``.

    :param str seq: The string to compare
    :param list rest: The list of strings to compare to
    :param int min_similarity: Minimum fraction to be considered similar

    :returns: If ``seq`` is similar to every sequence in ``rest``
    :rtype: bool

    """
    for comp_seq in rest:
        dist = dnautils.hamming(
            getattr(comp_seq, 'cdr3_' + field).replace('X', '-'),
            getattr(seq, 'cdr3_' + field).replace('X', '-')
        )
        sim_frac = 1 - dist / len(comp_seq.cdr3_aa)
        if sim_frac < min_similarity:
            return False
    return True


def can_subclone(sub_seqs, parent_seqs, min_similarity):
    for seq in sub_seqs:
        if not similar_to_all(seq, parent_seqs, min_similarity):
            return False
    return True


class ClonalWorker(concurrent.Worker):
    defaults = {
        # common
        'include_indels': False,
        'exclude_partials': False,
        'min_identity': 0.0,
        'min_copy': 2,
        'max_padding': None,

        # similarity
        'min_similarity': 0.85,

        # lineage
        'mut_cutoff': 4,
        'min_mut_occurrence': 2,
        'min_mut_samples': 1,
        'min_seq_instances': 1,
    }

    def __init__(self, session, **kwargs):
        self.session = session
        for prop, default in self.defaults.items():
            setattr(self, prop, kwargs.get(prop, default))
        for prop, value in kwargs.items():
            if prop not in self.defaults:
                setattr(self, prop, value)
        self._tasks = 0

    def get_bucket_seqs(self, bucket, sort):
        query = self.session.query(
            Sequence
        ).join(SequenceCollapse).filter(
            Sequence.subject_id == bucket.subject_id,
            Sequence.v_gene == bucket.v_gene,
            Sequence.j_gene == bucket.j_gene,
            Sequence.cdr3_num_nts == bucket.cdr3_num_nts,
            Sequence._insertions == bucket._insertions,
            Sequence._deletions == bucket._deletions,

            SequenceCollapse.copy_number_in_subject >= self.min_copy,
        )
        if sort:
            query = query.order_by(
                desc(SequenceCollapse.copy_number_in_subject),
                Sequence.ai
            )
        if self.min_identity > 0:
            query = query.filter(
                Sequence.v_match / Sequence.v_length >= self.min_identity
            )
        if self.min_seq_instances > 1:
            query = query.filter(
                SequenceCollapse.instances_in_subject >= self.min_seq_instances
            )
        if not self.include_indels:
            query = query.filter(Sequence.probable_indel_or_misalign == 0)
        if self.exclude_partials:
            query = query.filter(Sequence.partial == 0)
        if self.max_padding is not None:
            query = query.filter(Sequence.seq_start <= self.max_padding)
        return query

    def do_task(self, bucket):
        self.run_bucket(bucket)
        self._tasks += 1
        if self._tasks % 100 == 0:
            self.session.commit()
            self.info('Collapsed {} buckets'.format(self._tasks))

    def cleanup(self):
        self.session.commit()
        self.session.close()


class LineageClonalWorker(ClonalWorker):
    def run_bucket(self, bucket):
        updates = []
        consensus_needed = set([])

        seqs = self.get_bucket_seqs(bucket, sort=False)
        if seqs.count() > 0:
            cdr3_start = CDR3_OFFSET
            if bucket._insertions:
                cdr3_start += sum(
                    (i[1] for i in bucket._insertions)
                )
            germline = seqs[0].germline
            germline = ''.join((
                germline[:cdr3_start],
                '-' * bucket.cdr3_num_nts,
                germline[cdr3_start + bucket.cdr3_num_nts:]
            ))
            lineage = LineageWorker(
                self.session,
                clearcut.get_newick,
                self.min_mut_copies,
                self.min_mut_samples,
                exclude_stops=False,
                post_tree_hook=clearcut.minimize_tree
            )
            tree = lineage.get_tree(germline, seqs)

            for subtree in cut_tree(tree, self.mut_cutoff):
                new_clone = Clone(
                      subject_id=bucket.subject_id,
                      v_gene=bucket.v_gene,
                      j_gene=bucket.j_gene,
                      cdr3_num_nts=bucket.cdr3_num_nts,
                      _insertions=bucket._insertions,
                      _deletions=bucket._deletions
                )
                self.session.add(new_clone)
                self.session.flush()
                consensus_needed.add(new_clone.id)
                updates.extend([{
                    'sample_id': s[0],
                    'ai': s[1],
                    'clone_id': new_clone.id
                } for s in get_seq_pks(subtree)])

        if len(updates) > 0:
            self.session.bulk_update_mappings(Sequence, updates)
        generate_consensus(self.session, consensus_needed)


class SimilarityClonalWorker(ClonalWorker):
    def run_bucket(self, bucket):
        clones = OrderedDict()
        consensus_needed = set([])
        query = self.get_bucket_seqs(bucket, sort=True)

        if query.count() > 0:
            for seq in query:
                if seq.clone_id not in clones:
                    clones[seq.clone_id] = []
                clones[seq.clone_id].append(seq)
            if None in clones:
                for seq_to_add in clones[None]:
                    for clone_id, existing_seqs in clones.items():
                        if clone_id is None:
                            continue
                        if similar_to_all(seq_to_add, existing_seqs,
                                          self.level, self.min_similarity):
                            existing_seqs.append(seq_to_add)
                            break
                    else:
                        new_clone = Clone(subject_id=seq.subject_id,
                                          v_gene=seq.v_gene,
                                          j_gene=seq.j_gene,
                                          cdr3_num_nts=seq.cdr3_num_nts,
                                          _insertions=seq._insertions,
                                          _deletions=seq._deletions)
                        self.session.add(new_clone)
                        self.session.flush()
                        clones[new_clone.id] = [seq_to_add]
                del clones[None]

            for clone_id, seqs in clones.items():
                to_update = [
                    {
                        'sample_id': s.sample_id,
                        'ai': s.ai,
                        'clone_id': clone_id
                    } for s in seqs if s.clone_id is None
                ]
                if len(to_update) > 0:
                    self.session.bulk_update_mappings(Sequence, to_update)
                    consensus_needed.add(clone_id)
        generate_consensus(self.session, consensus_needed)


class SubcloneWorker(concurrent.Worker):
    def __init__(self, session, min_similarity):
        self.session = session
        self.min_similarity = min_similarity

    def do_task(self, bucket):
        clones = self.session.query(Clone).filter(
            Clone.subject_id == bucket.subject_id,
            Clone.v_gene == bucket.v_gene,
            Clone.j_gene == bucket.j_gene,
            Clone.cdr3_num_nts == bucket.cdr3_num_nts,
        ).all()

        if len(clones) == 0:
            return
        # The clones with indels are the only ones which can be subclones
        parent_clones = set([c for c in clones if len(c.insertions) == 0 and
                             len(c.deletions) == 0])
        potential_subclones = [c for c in clones if c not in parent_clones and
                               c.parent_id is None]
        self.info('Bucket {} has {} clones; parents={}, subs={}'.format(
            bucket, len(clones), len(parent_clones), len(potential_subclones)))

        for subclone in potential_subclones:
            for parent in parent_clones:
                if can_subclone(subclone.sequences, parent.sequences,
                                self.min_similarity):
                    subclone.parent = parent
                    break
        self.session.commit()

    def cleanup(self):
        self.session.commit()
        self.session.close()


def run_subclones(session, subject_ids, args):
    tasks = concurrent.TaskQueue()
    for subject_id in subject_ids:
        logger.info('Generating subclone task queue for subject {}'.format(
            subject_id))
        buckets = session.query(
            Clone.subject_id, Clone.v_gene, Clone.j_gene, Clone.cdr3_num_nts
        ).filter(
            Clone.subject_id == subject_id
        ).group_by(
            Clone.subject_id, Clone.v_gene, Clone.j_gene, Clone.cdr3_num_nts
        )
        for bucket in buckets:
            tasks.add_task(bucket)

    logger.info('Generated {} total subclone tasks'.format(tasks.num_tasks()))
    for i in range(0, min(tasks.num_tasks(), args.nproc)):
        tasks.add_worker(SubcloneWorker(config.init_db(args.db_config),
                                        args.similarity))
    tasks.start()


def run_clones(session, args):
    """Runs the clone-assignment pipeline stage.

    :param Session session: The database session
    :param Namespace args: The arguments passed to the command

    """
    if args.subject_ids is None:
        subject_ids = [s.id for s in session.query(Subject.id)]
    else:
        subject_ids = args.subject_ids
    mod_log.make_mod('clones', session=session, commit=True, info=vars(args))

    if args.regen:
        logger.info('Deleting existing clones')
        session.query(Clone).filter(
            Clone.subject_id.in_(subject_ids)
        ).delete(synchronize_session=False)
        session.commit()

    tasks = concurrent.TaskQueue()
    for subject_id in subject_ids:
        logger.info('Generating task queue for subject {}'.format(
            subject_id))
        buckets = session.query(
            Sequence.subject_id, Sequence.v_gene, Sequence.j_gene,
            Sequence.cdr3_num_nts, Sequence._insertions,
            Sequence._deletions
        ).filter(
            Sequence.subject_id == subject_id,
            Sequence.clone_id.is_(None)
        ).group_by(
            Sequence.subject_id, Sequence.v_gene, Sequence.j_gene,
            Sequence.cdr3_num_nts, Sequence._insertions,
            Sequence._deletions
        )
        for bucket in buckets:
            if not args.gene or bucket.v_gene.startswith(args.gene):
                tasks.add_task(bucket)

    logger.info('Generated {} total tasks'.format(tasks.num_tasks()))

    methods = {
        'similarity': SimilarityClonalWorker,
        'lineage': LineageClonalWorker,
    }
    for i in range(0, min(tasks.num_tasks(), args.nproc)):
        worker = methods[args.method](
            config.init_db(args.db_config), **args.__dict__
        )
        tasks.add_worker(worker)
    tasks.start()

    if args.subclones:
        run_subclones(session, subject_ids, args)
    else:
        logger.info('Skipping subclones')

    push_clone_ids(session)
    session.commit()
