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


def _consensus(strings):
    """Gets the unweighted consensus from a list of strings

    :param list strings: A set of equal-length strings.

    :returns: A consensus string
    :rtype: str

    """
    chrs = [Counter(chars).most_common(1)[0][0] for chars in zip(*strings)]
    return ''.join(chrs)


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

    def do_task(self, bucket):
        clones = OrderedDict()
        consensus_needed = set([])
        query = self._session.query(
            Sequence.sample_id, Sequence.ai, Sequence.clone_id,
            Sequence.cdr3_aa, Sequence.subject_id, Sequence.v_gene,
            Sequence.j_gene, Sequence.cdr3_num_nts
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
        query = query.order_by(
            desc(SequenceCollapse.copy_number_in_subject),
            Sequence.ai
        )

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
                        if self._similar_to_all(seq_to_add, existing_seqs):
                            existing_seqs.append(seq_to_add)
                            break
                    else:
                        new_clone = Clone(subject_id=seq.subject_id,
                                          v_gene=seq.v_gene,
                                          j_gene=seq.j_gene,
                                          cdr3_num_nts=seq.cdr3_num_nts)
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
        self._generate_consensus(consensus_needed)
        self._tasks += 1
        if self._tasks % 100 == 0:
            self._session.commit()
            self._print('Collapsed {} buckets'.format(self._tasks))

    def cleanup(self):
        self._session.commit()
        self._session.close()

    def _similar_to_all(self, seq, rest):
        """Determines if the string ``seq`` is at least ``min_similarity``
        similar to the list of strings ``rest``.

        :param str seq: The string to compare
        :param list rest: The list of strings to compare to

        :returns: If ``seq`` is similar to every sequence in ``rest``
        :rtype: bool

        """
        for comp_seq in rest:
            dist = dnautils.hamming(
                comp_seq.cdr3_aa.replace('X', '-'),
                seq.cdr3_aa.replace('X', '-')
            )
            sim_frac = 1 - dist / float(len(comp_seq.cdr3_aa))
            if sim_frac < self._min_similarity:
                return False
        return True

    def _generate_germline(self, seqs, clone):
        insertions = set([])
        for seq in seqs:
            if seq.insertions is not None:
                insertions.update(set(map(tuple, seq.insertions)))
        clone.insertions = insertions

        for seq in seqs:
            seq.clone_insertions = insertions

        rep_seq = seqs[0]
        rep_ins = rep_seq.insertions or 0
        if rep_ins != 0:
            rep_ins = sum((e[1] for e in rep_ins))
        germline = rep_seq.germline[:CDR3_OFFSET + rep_ins]

        for ins in insertions:
            if ins not in map(tuple, rep_seq.insertions):
                pos, size = ins
                germline = germline[:pos] + ('-' * size) + germline[pos:]
        germline += '-' * clone.cdr3_num_nts

        clone.functional = (
            len(germline) % 3 == 0 and
            not lookups.has_stop(germline)
        )

        j_region = rep_seq.germline.replace(
            '-', '')[-rep_seq.post_cdr3_length:]
        germline += j_region

        return germline

    def _generate_consensus(self, to_update):
        """Generates consensus CDR3s for clones.

        :param Session session: The database session
        :param list to_update: The list of clone IDs to assign to groups

        """

        if len(to_update) == 0:
            return
        for clone in funcs.periodic_commit(
                self._session,
                self._session.query(Clone).filter(Clone.id.in_(to_update)),
                interval=1000):
            seqs = self._session.query(
                Sequence
            ).join(SequenceCollapse).filter(
                Sequence.clone_id == clone.id,
                SequenceCollapse.copy_number_in_subject > 0
            ).all()
            clone.cdr3_nt = _consensus(map(lambda s: s.cdr3_nt, seqs))
            clone.cdr3_aa = lookups.aas_from_nts(clone.cdr3_nt)

            clone.germline = self._generate_germline(seqs, clone)

        self._session.commit()


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
        print 'Deleting existing clones'
        session.query(Clone).filter(
            Clone.subject_id.in_(subject_ids)
        ).delete(synchronize_session=False)
        session.commit()

    tasks = concurrent.TaskQueue()
    for subject_id in subject_ids:
        print 'Generating task queue for subject {}'.format(subject_id)
        buckets = session.query(
            Sequence.subject_id, Sequence.v_gene, Sequence.j_gene,
            Sequence.cdr3_num_nts, Sequence._insertions, Sequence._deletions
        ).filter(
            Sequence.subject_id == subject_id,
            Sequence.clone_id.is_(None)
        ).group_by(
            Sequence.subject_id, Sequence.v_gene, Sequence.j_gene,
            Sequence.cdr3_num_nts, Sequence._insertions, Sequence._deletions
        )
        for bucket in buckets:
            tasks.add_task(bucket)

    print 'Generated {} total tasks'.format(tasks.num_tasks())

    for i in range(0, min(tasks.num_tasks(), args.nproc)):
        tasks.add_worker(ClonalWorker(
            config.init_db(args.db_config),
            args.similarity / 100.0, args.include_indels,
            args.exclude_partials, args.min_identity / 100.0, args.min_copy,
            args.max_padding))
    tasks.start()

    print 'Pushing clone IDs to sample sequences'
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

    session.commit()
