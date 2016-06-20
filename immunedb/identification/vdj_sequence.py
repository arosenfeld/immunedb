import re

from Bio.Seq import Seq

import dnautils
from immunedb.common.models import CDR3_OFFSET
from immunedb.identification import AlignmentException, get_common_seq
from immunedb.identification.v_genes import VGene, find_v_position
from immunedb.util.funcs import find_streak_position
import immunedb.util.lookups as lookups


def gap_positions(seq):
    gaps = []
    for diff in re.finditer('[.]+', seq):
        start, end = diff.span()
        gaps.append((start, end - start))
    return gaps


class VDJSequence(object):
    MISMATCH_THRESHOLD = 3
    INDEL_WINDOW = 30
    INDEL_MISMATCH_THRESHOLD = .6

    __slots__ = [
        'ids', 'sequence', 'orig_sequence', 'v_germlines', 'j_germlines',
        '_force_vs', '_force_js', 'quality', 'orig_quality', '_j', '_j_start',
        'j_anchor_pos', 'j_anchor_len', 'j_length', 'j_match', '_v',
        'v_length', 'v_match', 'mutation_fraction', 'germline',
        'removed_prefix', 'removed_prefix_qual', '_cdr3_len', '_pad_len',
        'pre_cdr3_length', 'pre_cdr3_match', 'post_cdr3_length',
        'post_cdr3_match'
    ]

    def __init__(self, ids, seq, v_germlines, j_germlines,
                 force_vs=None, force_js=None, quality=None, analyze=False):
        self.ids = [ids] if type(ids) == str else ids
        self.sequence = seq.upper()
        self.orig_sequence = seq.upper()
        self.v_germlines = v_germlines
        self.j_germlines = j_germlines

        self._force_vs = force_vs
        self._force_js = force_js
        self.quality = quality
        self.orig_quality = quality

        self._j = None
        self.j_anchor_pos = None
        self.j_match = None

        self._v = None
        self.v_match = None

        self.mutation_fraction = None
        self.germline = None

        self.removed_prefix = ''
        self.removed_prefix_qual = ''
        if analyze:
            self.analyze()

    def analyze(self):
        if not all(map(lambda c: c in 'ATCGN', self.sequence)):
            raise AlignmentException('Invalid characters in sequence.')

        self._find_j()
        self._find_v()

    @property
    def j_gene(self):
        return self._j if self._j is None else sorted(self._j)

    @property
    def v_gene(self):
        return self._v if self._v is None else sorted(self._v)

    @v_gene.setter
    def v_gene(self, v):
        self._v = v

    @property
    def cdr3_start(self):
        return CDR3_OFFSET

    @property
    def num_gaps(self):
        return self.sequence[:self.cdr3_start].count('-')

    @property
    def cdr3(self):
        return self.sequence[self.cdr3_start:self.cdr3_start + self._cdr3_len]

    @property
    def partial(self):
        return self._pad_len > 0

    @property
    def pad_length(self):
        return self._pad_len if self._pad_len >= 0 else 0

    @property
    def in_frame(self):
        return len(self.cdr3) % 3 == 0 and self.cdr3_start % 3 == 0

    @property
    def stop(self):
        return lookups.has_stop(self.sequence)

    @property
    def functional(self):
        return self.in_frame and not self.stop

    def _check_j_with_missing(self, sequence, match):
        for pos in range(len(sequence) - len(match)):
            ss = sequence[pos:pos + len(match)]
            if dnautils.equal(ss, match):
                return pos
        return -1

    def _find_j(self):
        '''Finds the location and type of J gene'''
        # Iterate over every possible J anchor.  For each germline, try its
        # full sequence, then exclude the final 3 characters at a time until
        # there are only MIN_J_ANCHOR_LEN nucleotides remaining.
        #
        # For example, the order for one germline:
        # TGGTCACCGTCTCCTCAG
        # TGGTCACCGTCTCCT
        # TGGTCACCGTCT

        for match, j_gene in self.j_germlines.get_all_anchors(self._force_js):
            i = self.sequence.rfind(match)
            if i >= 0:
                return self._found_j(i, j_gene, match)

            rc = str(Seq(self.sequence).reverse_complement())
            i = rc.rfind(match)
            if i >= 0:
                self.sequence = rc
                if self.quality is not None:
                    self.quality = self.quality[::-1]
                return self._found_j(i, j_gene, match)

            i = self._check_j_with_missing(self.sequence, match)
            if i >= 0:
                return self._found_j(i, j_gene, match)

            i = self._check_j_with_missing(rc, match)
            if i >= 0:
                self.sequence = rc
                if self.quality is not None:
                    self.quality = self.quality[::-1]
                return self._found_j(i, j_gene, match)
        raise AlignmentException('Could not find J anchor')

    def _found_j(self, i, j_gene, match):
        # If a match is found, record its location and gene
        self.j_anchor_pos = i
        self.j_anchor_len = len(match)
        end_of_j = min(
            self.j_anchor_pos + self.j_germlines.anchor_len,
            len(self.sequence)
        )
        best_dist = None
        self._j = []
        if self._force_js:
            j_germs = {
                k: v for k, v in self.j_germlines.iteritems()
                if k in self._force_js
            }
        else:
            j_germs = self.j_germlines
        for j_gene, j_seq in j_germs.iteritems():
            seq_j = self.sequence[end_of_j - len(j_seq):end_of_j]
            dist = dnautils.hamming(seq_j, j_seq[:len(seq_j)])
            if best_dist is None or dist < best_dist:
                best_dist = dist
                self._j = set([j_gene])
            elif dist == best_dist:
                self._j.add(j_gene)

        if self._j is None:
            raise AlignmentException('Could not find suitable J anchor')

        # Get the full germline J gene
        j_full = self.j_germlines[self.j_gene[0]]

        # Get the portion of the germline J in the CDR3
        germline_in_cdr3 = self.j_germlines.get_j_in_cdr3(self.j_gene[0])
        cdr3_end_pos = (
            self.j_anchor_pos + self.j_germlines.anchor_len -
            self.j_germlines.upstream_of_cdr3
        )
        sequence_in_cdr3 = self.sequence[cdr3_end_pos - len(germline_in_cdr3):
                                         cdr3_end_pos]
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
        self._j_start = self.j_anchor_pos + len(match) - len(j_full)

        # If the trimmed germline J extends past the end of the
        # sequence, there is a misalignment
        if len(j_full) != len(
                self.sequence[self._j_start:self._j_start+len(j_full)]):
            self._j = None
            self.j_anchor_pos = None
            raise AlignmentException('Germline extended past end of J')

        self.j_length = len(j_full)

    def _find_v(self):
        for anchor_pos in find_v_position(self.sequence):
            self._found_v(anchor_pos)
            if self._v is not None:
                break

        if self._v is None:
            raise AlignmentException('Could not find suitable V anchor')

    def _found_v(self, anchor_pos):
        '''Finds the V gene closest to that of the sequence'''
        aligned_v = VGene(self.sequence)
        v_score = None
        for v, germ in sorted(self.v_germlines.alignments.iteritems()):
            if self._force_vs is not None and v not in self._force_vs:
                continue
            try:
                dist, total_length = germ.compare(aligned_v, self.j_anchor_pos,
                                                  self.MISMATCH_THRESHOLD)
            except:
                continue
            # Record this germline if it is has the lowest distance
            if dist is not None:
                if v_score is None or dist < v_score:
                    self._v = [v]
                    self.v_length = total_length
                    germ_pos = germ.ungapped_anchor_pos
                    v_score = dist
                elif dist == v_score:
                    # Add the V-tie
                    self._v.append(v)

        if self._v is not None:
            # Determine the pad length
            self._pad_len = germ_pos - anchor_pos
            # Mutation ratio is the distance divided by the length of overlap
            self.mutation_fraction = v_score / float(self.v_length)

    def align_to_germline(self, avg_len=None, avg_mut=None, trim_to=None):
        if avg_len is not None and avg_mut is not None:
            self._v = self.v_germlines.get_ties(self.v_gene, avg_len, avg_mut)
            self._j = self.j_germlines.get_ties(self.j_gene, avg_len, avg_mut)
        # Set the germline to the V gene up to the CDR3
        self.germline = get_common_seq(
            [self.v_germlines[v] for v in self._v]
        )[:CDR3_OFFSET]
        # If we need to pad the sequence, do so, otherwise trim the sequence to
        # the germline length
        if self._pad_len >= 0:
            self.sequence = 'N' * self._pad_len + str(self.sequence)
            if self.quality is not None:
                self.quality = (' ' * self._pad_len) + self.quality
        else:
            self.removed_prefix = self.sequence[:-self._pad_len]
            self.sequence = str(self.sequence[-self._pad_len:])
            if self.quality is not None:
                self.removed_prefix_qual = self.quality[:-self._pad_len]
                self.quality = self.quality[-self._pad_len:]
        # Update the anchor positions after adding padding / trimming
        self.j_anchor_pos += self._pad_len

        # Add germline gaps to sequence before CDR3 and update anchor positions
        for i, c in enumerate(self.germline):
            if c == '-':
                self.sequence = self.sequence[:i] + '-' + self.sequence[i:]
                if self.quality is not None:
                    self.quality = self.quality[:i] + ' ' + self.quality[i:]
                self.j_anchor_pos += 1

        j_germ = get_common_seq(
            map(reversed, [self.j_germlines[j] for j in self.j_gene]))
        j_germ = ''.join(reversed(j_germ))
        # Calculate the length of the CDR3
        self._cdr3_len = (
            self.j_anchor_pos + self.j_germlines.anchor_len -
            self.j_germlines.upstream_of_cdr3 - self.cdr3_start
        )

        if self._cdr3_len < 3:
            raise AlignmentException('CDR3 has no AAs'.format(self._cdr3_len))

        self.j_anchor_pos += self._cdr3_len
        # Fill germline CDR3 with gaps
        self.germline += '-' * self._cdr3_len
        self.germline += j_germ[-self.j_germlines.upstream_of_cdr3:]
        # If the sequence is longer than the germline, trim it
        if len(self.sequence) > len(self.germline):
            self.sequence = self.sequence[:len(self.germline)]
            if self.quality is not None:
                self.quality = self.quality[:len(self.germline)]
        elif len(self.sequence) < len(self.germline):
            self.sequence += 'N' * (len(self.germline) - len(self.sequence))
            if self.quality is not None:
                self.quality += ' ' * (len(self.germline) - len(self.quality))

        if trim_to is not None:
            old_padding = max(self._pad_len, 0)
            new_prefix = ''.join([
                c if c == '-' else 'N' for c in self.sequence[:trim_to]
            ])
            self.sequence = new_prefix + self.sequence[trim_to:]
            v_start = re.match('[N\-]*', self.sequence).span()[1]
            self._pad_len = self.sequence[:v_start].count('N')
            self.v_length -= self._pad_len - old_padding

        # Get the pre-CDR3 germline
        pre_cdr3_germ = self.germline[:self.cdr3_start]
        pre_cdr3_seq = self.sequence[:self.cdr3_start]

        # If there is padding, get rid of it in the sequence and align the
        # germline
        if self._pad_len > 0:
            pre_cdr3_germ = pre_cdr3_germ[self._pad_len:]
            pre_cdr3_seq = pre_cdr3_seq[self._pad_len:]

        # Calculate the pre-CDR3 length and distance
        self.pre_cdr3_length = len(pre_cdr3_seq)
        self.pre_cdr3_match = self.pre_cdr3_length - dnautils.hamming(
            str(pre_cdr3_seq), str(pre_cdr3_germ))

        # Get the length of J after the CDR3
        self.post_cdr3_length = self.j_germlines.upstream_of_cdr3
        # Get the sequence and germline sequences after CDR3
        post_j = j_germ[-self.post_cdr3_length:]
        post_s = self.sequence[-self.post_cdr3_length:]

        # Calculate their match count
        self.post_cdr3_match = self.post_cdr3_length - dnautils.hamming(
            post_j, post_s)

        self.v_match = self.v_length - dnautils.hamming(
            self.germline[:self.cdr3_start],
            self.sequence[:self.cdr3_start]
        )

        self.j_match = self.j_length - dnautils.hamming(
            self.germline[-len(j_germ):],
            self.sequence[-len(j_germ):]
        )

    @property
    def has_possible_indel(self):
        # Start comparison on first full AA to the INDEL_WINDOW or CDR3,
        # whichever comes first
        start = re.search('[ATCG]', self.sequence).start()
        germ = self.germline[start:self.cdr3_start]
        seq = self.sequence[start:self.cdr3_start]

        for i in range(0, len(germ) - self.INDEL_WINDOW + 1):
            dist = dnautils.hamming(germ[i:i+self.INDEL_WINDOW],
                                    seq[i:i+self.INDEL_WINDOW])
            if dist >= self.INDEL_MISMATCH_THRESHOLD * self.INDEL_WINDOW:
                return True

        return False
