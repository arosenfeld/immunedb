from collections import Counter, OrderedDict

from sqlalchemy import desc
from sqlalchemy.sql import text

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
        clone.cdr3_nt = consensus([s.cdr3_nt for s in seqs])
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


def consensus(strings):
    """Gets the unweighted consensus from a list of strings

    :param list strings: A set of equal-length strings.

    :returns: A consensus string
    :rtype: str

    """
    chrs = [Counter(chars).most_common(1)[0][0] for chars in zip(*strings)]
    return ''.join(chrs)


def similar_to_all(seq, rest, min_similarity):
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
            comp_seq.cdr3_aa.replace('X', '-'),
            seq.cdr3_aa.replace('X', '-')
        )
        sim_frac = 1 - dist / float(len(comp_seq.cdr3_aa))
        if sim_frac < min_similarity:
            return False
    return True


def can_subclone(sub_seqs, parent_seqs, min_similarity):
    for seq in sub_seqs:
        if not similar_to_all(seq, parent_seqs, min_similarity):
            return False
    return True


class ClonalWorker(concurrent.Worker):
    def __init__(self, session, min_similarity, include_indels,
                 exclude_partials, min_identity, min_copy, max_padding):
        self._session = session
        self._min_similarity = min_similarity
        self._include_indels = include_indels
        self._exclude_partials = exclude_partials
        self._min_identity = min_identity
        self._min_copy = min_copy
        self._max_padding = max_padding

        self._tasks = 0

    def get_query(self, bucket, sort):
        query = self._session.query(
            Sequence.sample_id, Sequence.ai, Sequence.clone_id,
            Sequence.cdr3_nt, Sequence.cdr3_aa, Sequence.subject_id,
            Sequence.v_gene, Sequence.j_gene, Sequence.cdr3_num_nts,
            Sequence._insertions, Sequence._deletions
        ).join(SequenceCollapse).filter(
            Sequence.subject_id == bucket.subject_id,
            Sequence.v_gene == bucket.v_gene,
            Sequence.j_gene == bucket.j_gene,
            Sequence.cdr3_num_nts == bucket.cdr3_num_nts,
            Sequence._insertions == bucket._insertions,
            Sequence._deletions == bucket._deletions,

            ~Sequence.cdr3_aa.like('%*%'),
            SequenceCollapse.copy_number_in_subject >= self._min_copy,
        )
        if sort:
            query = query.order_by(
                desc(SequenceCollapse.copy_number_in_subject),
                Sequence.ai
            )
        if self._min_identity > 0:
            query = query.filter(
                Sequence.v_match / Sequence.v_length >= self._min_identity
            )
        if not self._include_indels:
            query = query.filter(Sequence.probable_indel_or_misalign == 0)
        if self._exclude_partials:
            query = query.filter(Sequence.partial == 0)
        if self._max_padding is not None:
            query = query.filter(Sequence.pad_length <= self._max_padding)
        return query

    def do_task(self, args):
        if args['tcells']:
            self.tcell_clones(args['bucket'])
        else:
            self.bcell_clones(args['bucket'])

    def tcell_clones(self, bucket):
        updates = []
        clones = OrderedDict()
        consensus_needed = set([])

        for seq in self.get_query(bucket, False):
            key = (seq.v_gene, seq.j_gene, seq.cdr3_nt)
            if key in clones:
                clone = clones[key]
            else:
                for test_clone in clones.values():
                    same_bin = (test_clone.v_gene == key[0] and
                                test_clone.j_gene == key[1] and
                                test_clone.cdr3_num_nts == len(key[2]))
                    if same_bin and dnautils.equal(test_clone.cdr3_nt, key[2]):
                        clone = test_clone
                        break
                else:
                    new_clone = Clone(subject_id=seq.subject_id,
                                      v_gene=seq.v_gene,
                                      j_gene=seq.j_gene,
                                      cdr3_nt=seq.cdr3_nt,
                                      cdr3_num_nts=seq.cdr3_num_nts,
                                      _insertions=seq._insertions,
                                      _deletions=seq._deletions)
                    clones[key] = new_clone
                    self._session.add(new_clone)
                    self._session.flush()
                    clone = new_clone
                    consensus_needed.add(new_clone.id)
            updates.append({
                'sample_id': seq.sample_id,
                'ai': seq.ai,
                'clone_id': clone.id
            })

        if len(updates) > 0:
            self._session.bulk_update_mappings(Sequence, updates)
        generate_consensus(self._session, consensus_needed)

    def bcell_clones(self, bucket):
        clones = OrderedDict()
        consensus_needed = set([])
        query = self.get_query(bucket, True)

        if query.count() > 0:
            for seq in query:
                if seq.clone_id not in clones:
                    clones[seq.clone_id] = []
                clones[seq.clone_id].append(seq)
            if None in clones:
                for seq_to_add in clones[None]:
                    for clone_id, existing_seqs in clones.iteritems():
                        if clone_id is None:
                            continue
                        if similar_to_all(seq_to_add, existing_seqs,
                                          self._min_similarity):
                            existing_seqs.append(seq_to_add)
                            break
                    else:
                        new_clone = Clone(subject_id=seq.subject_id,
                                          v_gene=seq.v_gene,
                                          j_gene=seq.j_gene,
                                          cdr3_num_nts=seq.cdr3_num_nts,
                                          _insertions=seq._insertions,
                                          _deletions=seq._deletions)
                        self._session.add(new_clone)
                        self._session.flush()
                        clones[new_clone.id] = [seq_to_add]
                del clones[None]

            for clone_id, seqs in clones.iteritems():
                to_update = [
                    {
                        'sample_id': s.sample_id,
                        'ai': s.ai,
                        'clone_id': clone_id
                    } for s in seqs if s.clone_id is None
                ]
                if len(to_update) > 0:
                    self._session.bulk_update_mappings(Sequence, to_update)
                    consensus_needed.add(clone_id)
        generate_consensus(self._session, consensus_needed)
        self._tasks += 1
        if self._tasks % 100 == 0:
            self._session.commit()
            self.info('Collapsed {} buckets'.format(self._tasks))

    def cleanup(self):
        self._session.commit()
        self._session.close()


class SubcloneWorker(concurrent.Worker):
    def __init__(self, session, min_similarity):
        self._session = session
        self._min_similarity = min_similarity

    def do_task(self, bucket):
        clones = self._session.query(Clone).filter(
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
                                self._min_similarity):
                    subclone.parent = parent
                    break
        self._session.commit()

    def cleanup(self):
        self._session.commit()
        self._session.close()


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
                                        args.similarity / 100.0))
    tasks.start()


def run_clones(session, args):
    """Runs the clone-assignment pipeline stage.

    :param Session session: The database session
    :param Namespace args: The arguments passed to the command

    """
    if args.subject_ids is None:
        subject_ids = map(lambda s: s.id, session.query(Subject.id).all())
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
            tasks.add_task({
                'tcells': args.tcells,
                'bucket': bucket
            })

    logger.info('Generated {} total tasks'.format(tasks.num_tasks()))

    for i in range(0, min(tasks.num_tasks(), args.nproc)):
        tasks.add_worker(ClonalWorker(
            config.init_db(args.db_config),
            args.similarity / 100.0, args.include_indels,
            args.exclude_partials, args.min_identity / 100.0,
            args.min_copy, args.max_padding))
    tasks.start()

    if args.subclones:
        run_subclones(session, subject_ids, args)
    else:
        logger.info('Skipping subclones')

    push_clone_ids(session)
    session.commit()
