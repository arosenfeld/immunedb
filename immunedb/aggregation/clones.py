from collections import OrderedDict

from sqlalchemy import desc
from sqlalchemy.sql import text

import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import fcluster, linkage

import dnautils
from immunedb.common.models import (CDR3_OFFSET, Clone, Sequence,
                                    SequenceCollapse, Subject)
import immunedb.common.modification_log as mod_log
import immunedb.common.config as config
import immunedb.util.concurrent as concurrent
import immunedb.util.funcs as funcs
import immunedb.util.lookups as lookups
from immunedb.util.log import logger


def generate_consensus(session, clone_ids):
    """Generates consensus CDR3s for clones.

    :param Session session: The database session
    :param list clone_ids: The list of clone IDs to assign to groups

    """

    if not clone_ids:
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

        clone.functional = (
            clone.cdr3_num_nts % 3 == 0 and
            not lookups.has_stop(clone.germline)
        )

    session.commit()


def generate_germline(session, seqs, clone):
    rep_seq = seqs[0]
    germ = rep_seq.alignment_without_insertions[0]
    return ''.join((
        germ[:CDR3_OFFSET],
        '-' * clone.cdr3_num_nts,
        germ[CDR3_OFFSET + clone.cdr3_num_nts:]
    ))


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


def collapse_identical(session, buckets):
    for i, bucket in enumerate(buckets):
        clones = session.query(
            Clone.id, Clone.cdr3_aa
        ).filter(
            Clone.subject_id == bucket.subject_id,
            Clone.v_gene == bucket.v_gene,
            Clone.j_gene == bucket.j_gene,
            Clone.cdr3_num_nts == bucket.cdr3_num_nts,
        )
        if clones.count() < 2:
            continue
        logger.info('Reducing bucket {} / {} ({} clones)'.format(
            i, len(buckets), clones.count()))
        uniques = {}
        for c in clones:
            uniques.setdefault(c.cdr3_aa, []).append(c.id)
        uniques = [sorted(u) for u in uniques.values() if len(u) > 1]
        if len(uniques) > 1:
            logger.info('Collapsing {} duplicate CDR3s'.format(len(uniques)))
        for identical in uniques:
            rep_id = identical[0]
            session.query(Sequence).filter(
                Sequence.clone_id.in_(identical)
            ).update({'clone_id': rep_id}, synchronize_session=False)
            session.query(Clone).filter(
                Clone.id.in_(identical[1:])
            ).delete(synchronize_session=False)
    session.commit()


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


def can_subclone(sub_clone, parent_clone):
    return dnautils.equal(
        parent_clone.cdr3_aa.replace('X', '-'),
        sub_clone.cdr3_aa.replace('X', '-')
    )


class ClonalWorker(concurrent.Worker):
    defaults = {
        # common
        'exclude_partials': False,
        'min_copy': 2,
        'max_padding': None,
        'min_similarity': 0.85,
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

            SequenceCollapse.copy_number_in_subject >= self.min_copy,
        )
        if sort:
            query = query.order_by(
                desc(SequenceCollapse.copy_number_in_subject),
                Sequence.ai
            )
        if self.min_seq_instances > 1:
            query = query.filter(
                SequenceCollapse.instances_in_subject >= self.min_seq_instances
            )
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
                                          cdr3_num_nts=seq.cdr3_num_nts)
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


class ClusteringClonalWorker(ClonalWorker):
    def get_distances(self, seqs):
        dists = np.zeros((len(seqs), len(seqs)))
        for i, s1 in enumerate(seqs):
            for j, s2 in enumerate(seqs):
                d = dnautils.hamming(s1, s2)
                dists[i, j] = dists[j, i] = d / len(s1)
        return dists

    def assign_clones(self, df):
        seq = df.iloc[0]
        new_clone = Clone(subject_id=int(seq.subject_id),
                          v_gene=seq.v_gene,
                          j_gene=seq.j_gene,
                          cdr3_num_nts=int(seq.cdr3_num_nts))
        self.session.add(new_clone)
        self.session.flush()
        to_update = [
            {
                'sample_id': s.sample_id,
                'ai': s.ai,
                'clone_id': new_clone.id
            } for _, s in df.iterrows()
        ]
        self.session.bulk_update_mappings(Sequence, to_update)
        return int(new_clone.id)

    def run_bucket(self, bucket):
        seqs = self.get_bucket_seqs(bucket, sort=True)
        if seqs.count():
            df = pd.DataFrame({
                'sample_id': s.sample_id,
                'ai': s.ai,

                'subject_id': s.subject_id,
                'v_gene': s.v_gene,
                'j_gene': s.j_gene,
                'cdr3_num_nts': s.cdr3_num_nts,

                'cdr3': getattr(s, 'cdr3_' + self.level).replace('X', '-'),
                'current_clone_id': s.clone_id
            } for s in seqs)

            if len(df) == 1:
                df['clone'] = 1
            else:
                dists = squareform(self.get_distances(df.cdr3))
                clusters = fcluster(linkage(dists), self.min_similarity,
                                    criterion='distance')
                df['clone'] = clusters
            # NOTE: .groupby() doesn't work here since self.assign_clones has
            # side effects
            consensus_needed = []
            for clone in sorted(df.clone.unique()):
                new_id = self.assign_clones(df[df.clone == clone])
                consensus_needed.append(new_id)
            generate_consensus(self.session, consensus_needed)


class SubcloneWorker(concurrent.Worker):
    def __init__(self, session):
        self.session = session

    def do_task(self, bucket):
        clones = self.session.query(Clone).filter(
            Clone.subject_id == bucket.subject_id,
            Clone.v_gene == bucket.v_gene,
            Clone.j_gene == bucket.j_gene,
            Clone.cdr3_num_nts == bucket.cdr3_num_nts,
        ).all()

        if not clones:
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
                if can_subclone(subclone, parent):
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
        tasks.add_worker(SubcloneWorker(config.init_db(args.db_config)))
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

    if not args.skip_regen:
        logger.info('Deleting existing clones')
        session.query(Clone).filter(
            Clone.subject_id.in_(subject_ids)
        ).delete(synchronize_session=False)
        session.commit()

    tasks = concurrent.TaskQueue()
    all_buckets = []
    for subject_id in subject_ids:
        logger.info('Generating task queue for subject {}'.format(
            subject_id))
        buckets = session.query(
            Sequence.subject_id, Sequence.v_gene, Sequence.j_gene,
            Sequence.cdr3_num_nts
        ).filter(
            Sequence.subject_id == subject_id,
            Sequence.clone_id.is_(None)
        ).group_by(
            Sequence.subject_id, Sequence.v_gene, Sequence.j_gene,
            Sequence.cdr3_num_nts
        )
        for bucket in buckets:
            if not args.gene or bucket.v_gene.startswith(args.gene):
                tasks.add_task(bucket)
        all_buckets.extend(buckets)

    logger.info('Generated {} total tasks'.format(tasks.num_tasks()))

    methods = {
        'similarity': SimilarityClonalWorker,
        'cluster': ClusteringClonalWorker
    }
    for i in range(0, min(tasks.num_tasks(), args.nproc)):
        worker = methods[args.method](
            config.init_db(args.db_config), **args.__dict__
        )
        tasks.add_worker(worker)
    tasks.start()

    if args.skip_reduce:
        logger.info('Skipping reduce')
    else:
        collapse_identical(session, all_buckets)
    push_clone_ids(session)
    session.commit()

    if not args.skip_subclones:
        run_subclones(session, subject_ids, args)
    else:
        logger.info('Skipping subclones')
