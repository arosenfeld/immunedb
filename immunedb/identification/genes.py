import dnautils
import re

from Bio import SeqIO
from Bio.Seq import Seq

from immunedb.common.models import CDR3_OFFSET
from immunedb.util.hyper import hypergeom
from immunedb.identification import AlignmentException, get_common_seq
from immunedb.util.funcs import find_streak_position


class GermlineException(Exception):
    pass


class GeneName(object):
    def __init__(self, name):
        self.name = name
        parts = re.match('(([A-Z]+)(\d+)([^*]+)?)(\*(\d+))?',
                         self.name).groups()
        self.base = parts[0]
        self.prefix = parts[1]
        self.family = parts[2]
        self.allele = parts[5] if parts[5] else None

    def __str__(self):
        return self.name

    def __hash__(self):
        return hash(self.name)

    def __cmp__(self, other):
        return cmp(self.name, other.name)

    def __eq__(self, other):
        return hash(self) == hash(other)

    def __repr__(self):
        return ('<GeneName={}, base={}, prefix={}, family={}, '
                'allele={}>').format(str(self), self.base, self.prefix,
                                     self.family, self.allele)


class GeneTies(dict):
    TIES_PROB_THRESHOLD = 0.01

    def __init__(self, genes, remove_gaps=True):
        self.ties = {}
        self.hypers = {}
        self.remove_gaps = remove_gaps

        self.update(genes)

        self.allele_lookup = {}
        for name in self.keys():
            self.allele_lookup[name] = set([])
            for name2 in self.keys():
                if name2.base == name.base:
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
        length = int(length)
        mutation = round(mutation, 3)
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

            for name, v in sorted(self.iteritems()):
                s_2 = v.replace('-', '') if self.remove_gaps else v
                K = dnautils.hamming(s_1[-length:], s_2[-length:])
                p = self._hypergeom(length, mutation, K)
                if p >= self.TIES_PROB_THRESHOLD:
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


class VGermlines(GeneTies):
    def __init__(self, path_to_germlines):
        self._min_length = None
        self.alignments = {}

        with open(path_to_germlines) as fh:
            for record in SeqIO.parse(fh, 'fasta'):
                if record.seq.startswith('-'):
                    continue
                try:
                    v = VGene(str(record.seq))
                    name = GeneName(record.id)
                    self.alignments[name] = v
                    self[name] = str(record.seq)
                    if (self._min_length is None or
                            self._min_length > len(v.sequence_ungapped)):
                        self._min_length = len(v.sequence_ungapped)
                except Exception:
                    continue

        super(VGermlines, self).__init__({k: v for k, v in self.iteritems()})

    def get_single_tie(self, gene, length, mutation):
        return super(VGermlines, self).get_single_tie(
            gene, min(self.length_bucket(length), self._min_length), mutation
        )

    def length_bucket(self, length):
        if 0 < length <= 100:
            return 100
        if 100 < length <= 150:
            return 150
        if 150 < length <= 200:
            return 200
        return 300


class VGene(object):
    def __init__(self, gapped_sequence):
        self._gapped_seq = str(gapped_sequence).upper()
        if self._gapped_seq[CDR3_OFFSET:].count('-') > 0:
            raise AlignmentException('Cannot have gaps after CDR3 start '
                                     '(position {})'.format(CDR3_OFFSET))
        try:
            self._ungapped_anchor_pos = find_v_position(
                self.sequence_ungapped).next()
        except StopIteration:
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
        dist = dnautils.hamming(this_seq, other_seq)

        return dist, len(other_seq)


def find_v_position(sequence):
    '''Tries to find the end of the V gene region'''
    if type(sequence) == str:
        sequence = Seq(sequence)
    # Try to find DxxxyzC
    for found in _find_with_frameshifts(sequence, 'D(.{3}((YY)|(YC)|(YH)))C'):
        yield found
    # Try to find 'YYC', 'YCC', or 'YHC'
    for found in _find_with_frameshifts(sequence, 'Y([YHC])C'):
        yield found
    # Try to find 'DxxxxxC'
    for found in _find_with_frameshifts(sequence, 'D(.{5})C'):
        yield found


def _find_with_frameshifts(sequence, regex):
    for shift in [2, 1, 0]:
        seq = sequence[shift:]
        seq = seq[:len(seq) - len(seq) % 3]
        aas = str(seq.translate())
        res = re.search(regex, aas)
        if res is not None:
            yield (res.end() - 1) * 3 + shift


class JGermlines(GeneTies):
    defaults = {
        'upstream_of_cdr3': 31,
        'anchor_len': 18,
        'min_anchor_len': 12,
    }

    def __init__(self, path_to_germlines,
                 upstream_of_cdr3=defaults['upstream_of_cdr3'],
                 anchor_len=defaults['anchor_len'],
                 min_anchor_len=defaults['min_anchor_len']):
        self._upstream_of_cdr3 = upstream_of_cdr3
        self._anchor_len = anchor_len
        self._min_anchor_len = min_anchor_len
        self._min_length = None

        with open(path_to_germlines) as fh:
            for record in SeqIO.parse(fh, 'fasta'):
                name = GeneName(record.id)
                if all(map(lambda c: c in 'ATCGN', record.seq.upper())):
                    self[name] = str(record.seq).upper()
                    if (self._min_length is None or
                            len(self[name]) < self._min_length):
                        self._min_length = len(self[name])

        self._anchors = {name: seq[-anchor_len:] for name, seq in
                         self.iteritems()}
        super(JGermlines, self).__init__({k: v for k, v in self.iteritems()})

    @property
    def upstream_of_cdr3(self):
        return self._upstream_of_cdr3

    @property
    def anchor_len(self):
        return self._anchor_len

    @property
    def full_anchors(self):
        return self._anchors

    def get_j_in_cdr3(self, gene):
        return self[gene][:-self._upstream_of_cdr3]

    def get_all_anchors(self, allowed_genes=None):
        if allowed_genes is None:
            allowed_genes = self
        else:
            allowed_genes = {k: v for k, v in self.iteritems() if k.name in
                             allowed_genes}
        max_len = max(map(len, allowed_genes.values()))
        for trim_len in range(0, max_len, 3):
            for j, seq in allowed_genes.iteritems():
                trimmed_seq = seq[-self.anchor_len:-trim_len]
                if len(trimmed_seq) >= self._min_anchor_len:
                    yield trimmed_seq, j

    def get_single_tie(self, gene, length, mutation):
        seq = self[gene][-self.anchor_len:]
        tied = self.all_alleles(set([gene]))
        for j, other_seq in sorted(self.iteritems()):
            other_seq = other_seq[-self.anchor_len:][:len(seq)]
            if other_seq == seq:
                tied.add(j)
            elif dnautils.hamming(other_seq, seq) == 0:
                tied.add(j)
        return tied

    def all_ties(self, length, mutation):
        ties = {}
        for name in self:
            tie_name = tuple(sorted(self.get_ties([name], length, mutation)))
            if tie_name not in ties:
                ties[tie_name] = get_common_seq(
                    [self[n] for n in tie_name], right=True
                )
        return ties
