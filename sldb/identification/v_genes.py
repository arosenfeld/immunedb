import distance
import re
import itertools

import numpy as np
from scipy.stats import hypergeom

from Bio import SeqIO

import sldb.identification.anchors as anchors
from sldb.util.funcs import find_streak_position

def get_common_seq(seqs):
    v_gene = []
    for nts in itertools.izip_longest(*seqs, fillvalue='N'):
        v_gene.append(nts[0] if all(map(lambda n: n == nts[0], nts)) else 'N')
    v_gene = ''.join(v_gene)
    return v_gene[:VGene.CDR3_OFFSET]

class AlignmentException(Exception):
    pass

class VGene(object):
    CDR3_OFFSET = 309

    def __init__(self, gapped_sequence, force_anchor=False):
        self._gapped_seq = str(gapped_sequence)
        if self._gapped_seq[self.CDR3_OFFSET:].count('-') > 0:
            raise AlignmentException('Cannot have gaps after CDR3 start '
                                     '(position {})'.format(self.CDR3_OFFSET))
        if not force_anchor:
            self._ungapped_anchor_pos = anchors.find_v_position(
                self.sequence_ungapped)
        else:
            self._ungapped_anchor_pos = self.CDR3_OFFSET - \
                self._gapped_seq.count('-')
        if self._ungapped_anchor_pos is None:
            raise AlignmentException('Unable to find anchor')

    @property
    def sequence(self):
        return self._gapped_seq

    @property
    def sequence_ungapped(self):
        return self.sequence.replace('-', '')

    @property
    def ungapped_anchor_pos(self):
        return self._ungapped_anchor_pos

    def align(self, other_v):
        self._verify_type(other_v)

        diff = abs(self.ungapped_anchor_pos - other_v.ungapped_anchor_pos)
        this_seq = self.sequence_ungapped
        other_seq = other_v.sequence_ungapped

        # Trim the sequence which has the maximal anchor position, and
        # determine the CDR3 start position without gaps
        if self.ungapped_anchor_pos > other_v.ungapped_anchor_pos:
            this_seq = this_seq[diff:]
            cdr3_start = self.ungapped_anchor_pos - diff
        else:
            other_seq = other_seq[diff:]
            cdr3_start = other_v.ungapped_anchor_pos - diff

        return {
            'base': this_seq,
            'seq': other_seq,
            'diff': diff,
            'cdr3_start': cdr3_start
        }

    def compare(self, other_v, max_extent, max_streak):
        self._verify_type(other_v)

        alignment = self.align(other_v)
        this_seq = alignment['base'][:max_extent]
        other_seq = alignment['seq'][:max_extent]
        cdr3_offset = alignment['cdr3_start']

        # Determine the CDR3 in the germline and sequence
        this_cdr3 = this_seq[cdr3_offset:]
        other_cdr3 = other_seq[cdr3_offset:]
        length = min(len(this_cdr3), len(other_cdr3))
        this_cdr3 = this_cdr3[:length]
        other_cdr3 = other_cdr3[:length]
        if len(this_cdr3) == 0 or len(other_cdr3) == 0:
            raise AlignmentException('Empty CDR3 found after alignment')

        # Find the extent of the sequence's V into the CDR3
        streak = find_streak_position(
            this_cdr3, other_cdr3, max_streak)
        if streak is not None:
            # If there is a streak of mismatches, cut after the streak
            max_index = cdr3_offset + (streak - max_streak)
        else:
            # Unlikely: the CDR3 in the sequence exactly matches the
            # germline.  Use the smaller sequence length (full match)
            max_index = cdr3_offset + min(len(this_cdr3), len(other_cdr3))
        # Compare to the end of V
        this_seq = this_seq[:max_index]
        other_seq = other_seq[:max_index]

        if len(this_seq) != len(other_seq) or len(this_seq) == 0:
            raise AlignmentException('Unequal sequences after alignment')
        # Determine the distance between the germline and sequence
        dist = distance.hamming(this_seq, other_seq)

        return dist, len(other_seq)

    def _verify_type(self, other_v):
        if type(other_v) != type(self):
            raise AlignmentException('Must compare to instance of {}'.format(
                self.__class__.__name__))


class VGermlines(object):
    def __init__(self, path_to_germlines, ties_prob_threshold=.01):
        self._germlines = {}
        self._prob_threshold = ties_prob_threshold
        self._ties = {}
        self._min_length = None

        with open(path_to_germlines) as fh:
            for record in SeqIO.parse(fh, 'fasta'):
                if record.seq.startswith('-'):
                    continue
                try:
                    v = VGene(str(record.seq))
                    self._germlines[record.id] = v
                    if (self._min_length is None
                            or self._min_length > len(v.sequence_ungapped)):
                        self._min_length = len(v.sequence_ungapped)
                except:
                    pass

    def get_ties(self, genes, length, mutation):
        ties = set(genes)
        for gene in genes:
            ties.update(self._get_single_tie(gene, length, mutation))
        return ties

    def _get_single_tie(self, gene, length, mutation):
        length = min(self._min_length, self._length_bucket(length))
        mutation = self._mut_bucket(mutation)
        key = (length, mutation)

        if key not in self._ties:
            self._ties[key] = {}

        if gene not in self._ties[key]:
            s_1 = self._germlines[gene].sequence_ungapped
            self._ties[key][gene] = set([gene])

            for name, v in self._germlines.iteritems():
                s_2 = v.sequence_ungapped
                K = distance.hamming(s_1[-length:], s_2[-length:])
                dist = hypergeom(length, K, np.ceil(length * mutation))
                p = np.sum(
                    [dist.pmf(k) * np.power(.33, k) for k in xrange(int(np.ceil(K/2)), K)]
                )
                if p >= self._prob_threshold:
                    self._ties[key][gene].add(name)

        return self._ties[key][gene]

    def _length_bucket(self, length):
        if 0 < length <= 100:
            return 100
        if 100 < length <= 150:
            return 150
        if 150 < length <= 200:
            return 200
        return 300

    def _mut_bucket(self, mut):
        if 0 < mut <= .05:
            return .05
        if mut <= .15:
            return .15
        return .30

    def __getitem__(self, key):
        return self._germlines[key]

    def __iter__(self):
        for k in self._germlines.keys():
            yield k

    def iteritems(self):
        for k in self._germlines.keys():
            yield k, self._germlines[k]
