import re
import numpy as np

from Bio.Seq import Seq

import dnautils
from sldb.common.models import CDR3_OFFSET
from sldb.util.funcs import find_streak_position
from sldb.identification import AlignmentException
from sldb.identification.v_genes import VGene, get_common_seq, find_v_position
import sldb.identification.local_align as local_align


class VDJSequence(object):
    MISMATCH_THRESHOLD = 3
    INDEL_WINDOW = 30
    INDEL_MISMATCH_THRESHOLD = .6

    def __init__(self, id, seq, v_germlines, j_germlines,
                 force_vs=None, force_j=None, quality=None, skip_align=False):
        self._id = id
        self._seq = seq.upper()
        self.v_germlines = v_germlines
        self.j_germlines = j_germlines

        self._force_vs = force_vs
        self._force_j = force_j
        self._quality = quality

        self._j = None
        self._j_anchor_pos = None
        self._j_match = None

        self._v = None
        self._v_anchor_pos = None
        self._v_match = None

        self._mutation_frac = None
        self._germline = None
        self._cdr3_len = 0

        self._possible_indel = False

        if not all(map(lambda c: c in 'ATCGN', self.sequence)):
            raise AlignmentException('Invalid characters in sequence.')

        if not skip_align:
            self._find_j()
            self._find_v()

    def locally_align(self):
        vs, germ, new_seq = local_align.align_v(self._seq, self.v_germlines)

        self._v = set(vs)
        self._germline = germ[-CDR3_OFFSET:]
        self._seq = new_seq[len(germ) - CDR3_OFFSET:]
        if self._quality is not None:
            self._quality = self._quality[len(germ) - CDR3_OFFSET:]

            for i, c in enumerate(self._germline):
                if c == '-':
                    self._quality.insert(i, None)

        self._find_j(CDR3_OFFSET)
        self._cdr3_len = (
            self.j_anchor_pos + self.j_germlines.anchor_len -
             self.j_germlines.upstream_of_cdr3 - len(self._germline)
        )
        self._germline += ('-' * self._cdr3_len)

        j_germ = get_common_seq(
            map(reversed,
                [self.j_germlines[j][-self.j_germlines.upstream_of_cdr3:]
                    for j in self.j_gene]))
        self._germline += ''.join(
            reversed(j_germ[-self.j_germlines.upstream_of_cdr3:]))

        # TODO: These are not correct
        self._v_length = CDR3_OFFSET
        self._v_match = dnautils.hamming(self._germline[:CDR3_OFFSET],
                                         self._seq[:CDR3_OFFSET])
        self._pre_cdr3_length = self._v_length
        self._pre_cdr3_match = self._v_match
        self._post_cdr3_length = 0
        self._post_cdr3_match = 0

    @property
    def id(self):
        return self._id

    @property
    def quality(self):
        return self._quality

    @property
    def j_gene(self):
        if self._j is None:
            return None
        return sorted(self._j)

    @property
    def v_gene(self):
        if self._v is None:
            return None
        return sorted(self._v)

    @v_gene.setter
    def v_gene(self, v):
        self._v = v

    @property
    def j_anchor_pos(self):
        return self._j_anchor_pos

    @property
    def v_anchor_pos(self):
        return self._v_anchor_pos

    @property
    def cdr3(self):
        return self.sequence[CDR3_OFFSET:CDR3_OFFSET + self._cdr3_len]

    @property
    def partial(self):
        return self._pad_len > 0

    @property
    def sequence(self):
        return self._seq

    @property
    def aligned_v(self):
        return self._aligned_v

    @property
    def germline(self):
        return self._germline

    @property
    def mutation_fraction(self):
        return self._mutation_frac

    @property
    def num_gaps(self):
        return self.sequence[:CDR3_OFFSET].count('-')

    @property
    def pad_length(self):
        return self._pad_len if self._pad_len >= 0 else 0

    @property
    def in_frame(self):
        return len(self.cdr3) % 3 == 0

    @property
    def stop(self):
        for i in range(0, len(self.sequence), 3):
            if self.sequence[i:i+3] in ['TAG', 'TAA', 'TGA']:
                return True
        return False

    @property
    def functional(self):
        return self.in_frame and not self.stop

    @property
    def j_length(self):
        return self._j_length

    @property
    def j_match(self):
        return self._j_match

    @property
    def v_length(self):
        return self._v_length

    @property
    def v_match(self):
        return self._v_match

    @property
    def pre_cdr3_length(self):
        return self._pre_cdr3_length

    @property
    def pre_cdr3_match(self):
        return self._pre_cdr3_match

    @property
    def post_cdr3_length(self):
        return self._post_cdr3_length

    @property
    def post_cdr3_match(self):
        return self._post_cdr3_match

    def _find_j(self, offset=0):
        '''Finds the location and type of J gene'''
        # Iterate over every possible J anchor.  For each germline, try its
        # full sequence, then exclude the final 3 characters at a time until
        # there are only MIN_J_ANCHOR_LEN nucleotides remaining.
        #
        # For example, the order for one germline:
        # TGGTCACCGTCTCCTCAG
        # TGGTCACCGTCTCCT
        # TGGTCACCGTCT

        for match, full_anchor, j_gene in self.j_germlines.get_all_anchors(
                limit_genes=self._force_j):
            i = self._seq[offset:].rfind(match) + offset
            if i >= 0:
                return self._found_j(i, j_gene, match, full_anchor)

            i = Seq(self._seq).reverse_complement().rfind(match)
            if i >= 0:
                self._seq = str(Seq(self._seq).reverse_complement())
                if self._quality is not None:
                    self._quality.reverse()
                return self._found_j(i, j_gene, match, full_anchor)
        raise AlignmentException('Could not find J anchor')

    def _found_j(self, i, j_gene, match, full_anchor):
        # If a match is found, record its location and gene
        self._j_anchor_pos = i
        self._j = [j_gene]

        # Get the full germline J gene
        j_full = self.j_germlines[self.j_gene[0]]

        # Get the portion of the germline J in the CDR3
        germline_in_cdr3 = self.j_germlines.get_j_in_cdr3(self.j_gene[0])
        cdr3_end_pos = (
            self._j_anchor_pos + self.j_germlines.anchor_len -
            self.j_germlines.upstream_of_cdr3
        )
        sequence_in_cdr3 = self.sequence[cdr3_end_pos - len(germline_in_cdr3)
            :cdr3_end_pos]
        if len(germline_in_cdr3) == 0 or len(sequence_in_cdr3) == 0:
            self._j = None
            raise AlignmentException('Could not find sequence or germline in '
                                     'CDR3')

        # Get the extent of the J in the CDR3
        streak = find_streak_position(
            reversed(germline_in_cdr3),
            reversed(sequence_in_cdr3),
            self.MISMATCH_THRESHOLD)

        # Trim the J gene based on the extent in the CDR3
        if streak is not None:
            j_full = j_full[len(germline_in_cdr3) - streak:]

        # Find where the full J starts
        self._j_start = self._j_anchor_pos + len(match) - len(j_full)

        # If the trimmed germline J extends past the end of the
        # sequence, there is a misalignment
        if len(j_full) != len(
                self.sequence[self._j_start:self._j_start+len(j_full)]):
            self._j = None
            self._j_anchor_pos = None
            raise AlignmentException('Germline extended past end of J')

        # Get the full-J distance
        dist = dnautils.hamming(
            j_full,
            str(self.sequence[self._j_start:self._j_start+len(j_full)])
        )

        self._j = self.j_germlines.get_ties(self.j_gene[0], match)
        self._j_length = len(j_full)

    def _find_v(self):
        for self._v_anchor_pos in find_v_position(self.sequence):
            self._found_v()
            if self._v is not None:
                break

        if self._v is None:
            raise AlignmentException('Could not find suitable V anchor')

    def _found_v(self):
        '''Finds the V gene closest to that of the sequence'''
        self._v_score = None
        self._aligned_v = VGene(self.sequence)
        if self._force_vs is not None:
            germlines = {vg: self.v_germlines[vg] for vg in self._force_vs}
        else:
            germlines = self.v_germlines

        for v, germ in sorted(germlines.iteritems()):
            try:
                dist, total_length = germ.compare(self._aligned_v,
                                                  self._j_anchor_pos,
                                                  self.MISMATCH_THRESHOLD)
            except:
                continue
            # Record this germline if it is has the lowest distance
            if dist is not None:
                if self._v_score is None or dist < self._v_score:
                    self._v = [v]
                    self._v_score = dist
                    self._v_length = total_length
                    self._germ_pos = germ.ungapped_anchor_pos
                elif dist == self._v_score:
                    # Add the V-tie
                    self._v.append(v)

        if self._v is not None:
            # Determine the pad length
            self._pad_len = self._germ_pos - self.v_anchor_pos
            # Mutation ratio is the distance divided by the length of overlap
            self._mutation_frac = self._v_score / float(self._v_length)

    def align_to_germline(self, avg_len=None, avg_mut=None):
        if avg_len is not None and avg_mut is not None:
            self._v = self.v_germlines.get_ties(self.v_gene, avg_len, avg_mut)
        # Set the germline to the V gene up to the CDR3
        self._germline = get_common_seq(
            [self.v_germlines[v].sequence for v in self._v]
        )[:CDR3_OFFSET]
        self._pad_len = (len(self._germline.replace('-', '')) -
                         self.v_anchor_pos)
        # If we need to pad the sequence, do so, otherwise trim the sequence to
        # the germline length
        if self._pad_len >= 0:
            self._seq = 'N' * self._pad_len + str(self._seq)
            if self._quality is not None:
                self._quality = ([None] * self._pad_len) + self._quality
        else:
            self._seq = str(self._seq[-self._pad_len:])
            if self._quality is not None:
                self._quality = self._quality[-self._pad_len:]
        # Update the anchor positions after adding padding / trimming
        self._j_anchor_pos += self._pad_len
        self._v_anchor_pos += self._pad_len

        # Add germline gaps to sequence before CDR3 and update anchor positions
        for i, c in enumerate(self._germline):
            if c == '-':
                self._seq = self._seq[:i] + '-' + self._seq[i:]
                if self._quality is not None:
                    self._quality.insert(i, None)
                self._j_anchor_pos += 1
                self._v_anchor_pos += 1

        j_germ = get_common_seq(
            map(reversed, [self.j_germlines[j] for j in self.j_gene]))
        j_germ = ''.join(reversed(j_germ))
        # Calculate the length of the CDR3
        self._cdr3_len = (
            self.j_anchor_pos + self.j_germlines.anchor_len -
             self.j_germlines.upstream_of_cdr3 - self.v_anchor_pos
        )

        if self._cdr3_len <= 0:
            self._v = None
            raise AlignmentException('CDR3 has non-positive length: {}'
                                     '.'.format(self._cdr3_len))

        self._j_anchor_pos += self._cdr3_len
        # Fill germline CDR3 with gaps
        self._germline += '-' * self._cdr3_len
        self._germline += j_germ[-self.j_germlines.upstream_of_cdr3:]
        # If the sequence is longer than the germline, trim it
        if len(self.sequence) > len(self.germline):
            self._seq = self._seq[:len(self._germline)]
            if self._quality is not None:
                self._quality = self._quality[:len(self._germline)]
        elif len(self.sequence) < len(self.germline):
            # If the germline is longer than the sequence, there was probably a
            # deletion, so flag it as such
            self._seq += 'N' * (len(self.germline) - len(self.sequence))
            if self._quality is not None:
                self._quality.extend([None] * (len(self.germline) -
                                     len(self._quality)))
            self._possible_indel = True

        # Get the pre-CDR3 germline and sequence stripped of gaps
        pre_cdr3_germ = self.germline[:CDR3_OFFSET].replace('-', '')
        pre_cdr3_seq = self.sequence[:CDR3_OFFSET].replace('-', '')
        # If there is padding, get rid of it in the sequence and align the
        # germline
        if self._pad_len > 0:
            pre_cdr3_germ = pre_cdr3_germ[self._pad_len:]
            pre_cdr3_seq = pre_cdr3_seq[self._pad_len:]

        # Calculate the pre-CDR3 length and distance
        self._pre_cdr3_length = len(pre_cdr3_seq)
        self._pre_cdr3_match = self._pre_cdr3_length - dnautils.hamming(
            str(pre_cdr3_seq), str(pre_cdr3_germ))

        # Get the length of J after the CDR3
        self._post_cdr3_length = self.j_germlines.upstream_of_cdr3
        # Get the sequence and germline sequences after CDR3
        post_j = j_germ[-self.post_cdr3_length:]
        post_s = self.sequence[CDR3_OFFSET+len(self.cdr3):]

        # Calculate their match count
        self._post_cdr3_match = self.post_cdr3_length - dnautils.hamming(
            post_j, post_s)

        # Updaself.sequence[:CDR3_OFFSET]te the V and J matches after ties
        self._v_match = self.v_length - dnautils.hamming(
            self.germline[:CDR3_OFFSET],
            self.sequence[:CDR3_OFFSET]
        )

        self._j_match = self.j_length - dnautils.hamming(
            self.germline[-len(j_germ):],
            self.sequence[-len(j_germ):]
        )

    @property
    def has_possible_indel(self):
        # Start comparison on first full AA to the INDEL_WINDOW or CDR3,
        # whichever comes first
        if self._possible_indel:
            return True

        start = re.search('[ATCG]', self.sequence).start()
        germ = self.germline[start:CDR3_OFFSET]
        seq = self.sequence[start:CDR3_OFFSET]

        for i in range(0, len(germ) - self.INDEL_WINDOW + 1):
            dist = dnautils.hamming(germ[i:i+self.INDEL_WINDOW],
                                    seq[i:i+self.INDEL_WINDOW])
            if dist >= self.INDEL_MISMATCH_THRESHOLD * self.INDEL_WINDOW:
                return True

        return False
