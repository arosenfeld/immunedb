import dnautils
import itertools
import traceback

from immunedb.common.models import (CDR3_OFFSET, DuplicateSequence, NoResult,
                                    Sequence)
import immunedb.util.funcs as funcs
import immunedb.util.lookups as lookups
from immunedb.util.hyper import hypergeom
from immunedb.util.log import logger


class AlignmentException(Exception):
    pass


def add_as_noresult(session, vdj, sample, reason):
    try:
        session.bulk_save_objects([
            NoResult(
                seq_id=seq_id,
                sample_id=sample.id,
                sequence=vdj.orig_sequence,
                quality=vdj.orig_quality,
                reason=reason
            ) for seq_id in vdj.ids
        ])
    except ValueError:
        pass


def add_as_sequence(session, vdj, sample):
    try:
        seq = Sequence(
            seq_id=vdj.ids[0],
            sample_id=sample.id,

            subject_id=sample.subject.id,

            partial=vdj.partial,

            probable_indel_or_misalign=vdj.has_possible_indel,

            v_gene=funcs.format_ties(vdj.v_gene, vdj.v_germlines.prefix,
                                     strip_alleles=True),
            j_gene=funcs.format_ties(vdj.j_gene, vdj.j_germlines.prefix,
                                     strip_alleles=True),

            num_gaps=vdj.num_gaps,
            pad_length=vdj.pad_length,

            v_match=vdj.v_match,
            v_length=vdj.v_length,
            j_match=vdj.j_match,
            j_length=vdj.j_length,

            removed_prefix=vdj.removed_prefix,
            removed_prefix_qual=vdj.removed_prefix_qual,
            v_mutation_fraction=vdj.mutation_fraction,

            pre_cdr3_length=vdj.pre_cdr3_length,
            pre_cdr3_match=vdj.pre_cdr3_match,
            post_cdr3_length=vdj.post_cdr3_length,
            post_cdr3_match=vdj.post_cdr3_match,

            in_frame=vdj.in_frame,
            functional=vdj.functional,
            stop=vdj.stop,
            copy_number=len(vdj.ids),

            cdr3_nt=vdj.cdr3,
            cdr3_num_nts=len(vdj.cdr3),
            cdr3_aa=lookups.aas_from_nts(vdj.cdr3),

            sequence=str(vdj.sequence),
            quality=vdj.quality,

            germline=vdj.germline)
        session.add(seq)
        session.flush()

        # Add duplicate sequences
        try:
            session.bulk_save_objects([
                DuplicateSequence(
                    sample_id=sample.id,
                    seq_id=seq_id,
                    duplicate_seq_ai=seq.ai
                ) for seq_id in vdj.ids[1:]
            ])
        except ValueError:
            pass
    except ValueError as e:
        add_as_noresult(session, vdj, sample, str(e))


def add_uniques(session, sample, vdjs, realign_len=None,
                realign_mut=None, min_similarity=0, max_vties=50,
                trim_to=None, max_padding=None):
    bucketed_seqs = {}
    for vdj in funcs.periodic_commit(session, vdjs):
        try:
            if realign_len is not None:
                vdj.align_to_germline(realign_len, realign_mut, trim_to)
            if vdj.v_match / float(vdj.v_length) < min_similarity:
                raise AlignmentException(
                    'V-identity too low {} < {}'.format(
                        vdj.v_match / float(vdj.v_length), min_similarity))
            if len(vdj.v_gene) > max_vties:
                raise AlignmentException('Too many V-ties {} > {}'.format(
                    len(vdj.v_gene), max_vties))
            if max_padding is not None and vdj.pad_length > max_padding:
                raise AlignmentException('Too much padding {} (max {})'.format(
                    vdj.pad_length, max_padding
                ))
            bucket_key = (
                funcs.format_ties(vdj.v_gene, vdj.v_germlines.prefix,
                                  strip_alleles=True),
                funcs.format_ties(vdj.j_gene, vdj.j_germlines.prefix,
                                  strip_alleles=True),
                len(vdj.cdr3)
            )
            if bucket_key not in bucketed_seqs:
                bucketed_seqs[bucket_key] = {}
            bucket = bucketed_seqs[bucket_key]

            if vdj.sequence in bucket:
                bucket[vdj.sequence].ids += vdj.ids
            else:
                bucket[vdj.sequence] = vdj
        except AlignmentException as e:
            add_as_noresult(session, vdj, sample, str(e))
        except:
            logger.error('\tUnexpected error processing sequence '
                         '{}\n\t{}'.format(vdj.ids[0], traceback.format_exc()))

    # Collapse sequences that are the same except for Ns
    for sequences in funcs.periodic_commit(session, bucketed_seqs.values()):
        sequences = sorted(sequences.values(), cmp=lambda a, b:
                           cmp(len(a.ids), len(b.ids)))
        while len(sequences) > 0:
            larger = sequences.pop(0)
            for i in reversed(range(len(sequences))):
                smaller = sequences[i]

                if dnautils.equal(larger.sequence, smaller.sequence):
                    larger.ids += smaller.ids
                    del sequences[i]
            add_as_sequence(session, larger, sample)
    session.commit()


class GeneTies(dict):
    def __init__(self, genes, remove_gaps=True, ties_prob_threshold=.01):
        self.ties = {}
        self.hypers = {}
        self.ties_prob_threshold = ties_prob_threshold
        self.remove_gaps = remove_gaps

        self.update(genes)

        self.allele_lookup = {}
        for name in self.keys():
            self.allele_lookup[name] = set([])
            for name2 in self.keys():
                if name2.split('*')[0] == name.split('*')[0]:
                    self.allele_lookup[name].add(name2)

    def all_ties(self, length, mutation, cutoff=True):
        ties = {}
        for name in self:
            tie_name = tuple(sorted(self.get_ties([name], length, mutation)))
            if tie_name not in ties:
                ties[tie_name] = get_common_seq(
                    [self[n] for n in tie_name], cutoff=cutoff
                )
        return ties

    def get_ties(self, genes, length, mutation):
        ties = set([])
        for gene in genes:
            ties.update(self.get_single_tie(gene, length, mutation))
        return ties

    def get_single_tie(self, gene, length, mutation):
        mutation = self.mut_bucket(mutation)
        key = (length, mutation)

        if key not in self.ties:
            self.ties[key] = {}

        if gene not in self:
            return set([gene])

        if gene not in self.ties[key]:
            s_1 = (
                self[gene].replace('-', '') if self.remove_gaps else self[gene]
            )
            self.ties[key][gene] = set([gene])

            for name, v in self.iteritems():
                s_2 = v.replace('-', '') if self.remove_gaps else v
                K = dnautils.hamming(s_1[-length:], s_2[-length:])
                p = self._hypergeom(length, mutation, K)
                if p >= self.ties_prob_threshold:
                    self.ties[key][gene].add(name)
            self.ties[key][gene] = self.all_alleles(self.ties[key][gene])

        return self.ties[key][gene]

    def _hypergeom(self, length, mutation, K):
        key = (length, mutation, K)
        if key not in self.hypers:
            self.hypers[key] = hypergeom(length, mutation, K)
        return self.hypers[key]

    def mut_bucket(self, mut):
        if 0 <= mut <= .05:
            return .05
        if mut <= .15:
            return .15
        return .30

    def all_alleles(self, genes):
        all_genes = set([])
        for gene in genes:
            all_genes.update(self.allele_lookup[gene])
        return all_genes


def get_common_seq(seqs, cutoff=True):
    if len(seqs) == 0:
        return seqs[0]
    v_gene = []
    for nts in itertools.izip_longest(*seqs, fillvalue='N'):
        v_gene.append(nts[0] if all(map(lambda n: n == nts[0], nts)) else 'N')
    v_gene = ''.join(v_gene)
    if cutoff:
        return v_gene[:CDR3_OFFSET]
    return v_gene
