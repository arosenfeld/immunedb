import dnautils
import re

from Bio import SeqIO
from Bio.Seq import Seq

from immunedb.common.models import CDR3_OFFSET
from immunedb.identification import AlignmentException, GeneTies
from immunedb.util.funcs import find_streak_position


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
        dist = dnautils.hamming(this_seq, other_seq)

        return dist, len(other_seq)

    def _verify_type(self, other_v):
        if type(other_v) != type(self):
            raise AlignmentException('Must compare to instance of {}'.format(
                self.__class__.__name__))


class VGermlines(GeneTies):
    ALLOWED_PREFIX = set([
        'IGHV',
        'TRAV',
        'TRBV',
        'TRDV',
        'TRGV',
    ])
    def __init__(self, path_to_germlines, ties_prob_threshold=.01,
                 include_prepadded=False):
        self._min_length = None
        self.alignments = {}

        with open(path_to_germlines) as fh:
            for record in SeqIO.parse(fh, 'fasta'):
                if record.seq.startswith('-') and not include_prepadded:
                    continue
                self.prefix = record.id[:4]
                assert self.prefix in VGermlines.ALLOWED_PREFIX
                try:
                    v = VGene(str(record.seq))
                    self.alignments[record.id] = v
                    self[record.id] = str(record.seq)
                    if (self._min_length is None or
                            self._min_length > len(v.sequence_ungapped)):
                        self._min_length = len(v.sequence_ungapped)
                except:
                    continue

        super(VGermlines, self).__init__(
            {k: v for k, v in self.iteritems()},
            ties_prob_threshold=ties_prob_threshold
        )

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
