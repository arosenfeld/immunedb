import re
import numpy as np

from Bio.Seq import Seq

import dnautils
from sldb.common.models import CDR3_OFFSET
from sldb.util.funcs import find_streak_position
from sldb.identification import AlignmentException
from sldb.identification.v_genes import VGene, get_common_seq, find_v_position


def get_gap_differences(reference, seq):
    diffs = []
    for diff in re.finditer('[.]+', reference):
        start, end = diff.span()
        diffs.append((start, end - start))

    return diffs


class VDJSequence(object):
    MISMATCH_THRESHOLD = 3
    INDEL_WINDOW = 30
    INDEL_MISMATCH_THRESHOLD = .6

    def __init__(self, id, seq, v_germlines, j_germlines,
                 force_vs=None, force_j=None, quality=None,
                 locally_align=False):
        self.id = id
        self.sequence = seq.upper()
        self.v_germlines = v_germlines
        self.j_germlines = j_germlines

        self._force_vs = force_vs
        self._force_j = force_j
        self.quality = quality

        self._j = None
        self.j_anchor_pos = None
        self.j_match = None

        self._v = None
        self._v_end = None
        self.v_match = None

        self.mutation_fraction = None
        self.germline = None
        self._cdr3_len = 0

        self._possible_indel = False
        self.regions = [78, 36, 51, 30, 114]

        self.insertions = None
        self.deletions = None

        self.removed_prefix = ''
        self.removed_prefix_qual = ''

        if not all(map(lambda c: c in 'ATCGN', self.sequence)):
            raise AlignmentException('Invalid characters in sequence.')

        if locally_align:
            self._locally_align(*locally_align)
        else:
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
    def cdr3(self):
        return self.sequence[self._v_end:self._v_end + self._cdr3_len]

    @property
    def partial(self):
        return self._pad_len > 0

    @property
    def num_gaps(self):
        return self.sequence[:self._v_end].count('-')

    @property
    def pad_length(self):
        return self._pad_len if self._pad_len >= 0 else 0

    @property
    def in_frame(self):
        return len(self.cdr3) % 3 == 0 and self._v_end % 3 == 0

    @property
    def stop(self):
        for i in range(0, len(self.sequence), 3):
            if self.sequence[i:i+3] in ['TAG', 'TAA', 'TGA']:
                return True
        return False

    @property
    def functional(self):
        return self.in_frame and not self.stop

    def _locally_align(self, avg_mut, avg_len, insert_penalty=-50,
                      delete_penalty=-50, extend_penalty=-10,
                      mismatch_penalty=-20, match_score=40):
        max_align = None
        for name, germr in self.v_germlines.iteritems():
            germ = germr.sequence.replace('-', '')
            try:
                g, s, germ_omit, seq_omit, score = dnautils.align(
                    germ, self.sequence, insert_penalty, delete_penalty,
                    extend_penalty, mismatch_penalty, match_score)
                s = ''.join(reversed(s))
                g = ''.join(reversed(g))
            except Exception as e:
                continue
            if max_align is None or score >= max_align['score']:
                v_length = len(s)
                q = None
                if self.quality is not None:
                    q = self.quality[seq_omit:]
                if len(germ) > len(g):
                    g = germ[:germ_omit].lower() + g
                    g += germ[len(g) - g.count('-'):]

                if len(self.sequence) > len(s):
                    s += self.sequence[seq_omit + len(s) - s.count('-'):]

                max_align = {
                    'original_germ': germr.sequence,
                    'seq': s,
                    'qual': q,
                    'germ': g,
                    'score': score,
                    'v_length': v_length,
                    'vs': [(name, g)]
                }
            elif score == max_align['score']:
                max_align['vs'].append((name, g))

        if max_align is None:
            raise AlignmentException('Could not locally align sequence.')

        self._v = map(lambda e: e[0], max_align['vs'])

        imgt_gaps = [
            i for i, c in enumerate(max_align['original_germ']) if c == '-'
        ]

        self.germline = list(
            get_common_seq(
                map(lambda e: e[1], max_align['vs']), cutoff=False
            ).replace('-', '.')
        )
        self.sequence = list(max_align['seq'].replace('-', '.'))

        if self.quality is not None:
            self.quality = max_align['qual']
            for i, c in enumerate(self.sequence):
                if c == '.':
                    self.quality.insert(i, None)

        for gap in imgt_gaps:
            self.germline.insert(gap, '-')
            self.sequence.insert(gap, '-')
            if self.quality is not None:
                self.quality.insert(gap, None)

        self.germline = ''.join(self.germline)
        self.sequence = ''.join(self.sequence)

        offset = re.search('[ATCGN]', self.germline)
        if offset is None:
            raise AlignmentException('Entire germline gapped.')

        offset = offset.start()
        self.germline = self.germline[offset:]
        self.sequence = self.sequence[offset:]
        if self.quality is not None:
            self.quality = self.quality[offset:]

        self.v_length = max_align['v_length']

        # TODO: Add padding for partials
        self._pad_len = 0
        if self.v_length is None:
            raise AlignmentException('Germline was too small.')

        self.insertions = get_gap_differences(self.germline, self.sequence)
        self.deletions = get_gap_differences(self.sequence, self.germline)

        self._v = self.v_germlines.get_ties(self.v_gene, avg_len, avg_mut)
        insertion_total = sum(map(lambda e: e[1], self.insertions))
        common_v = get_common_seq(
            [self.v_germlines[v].sequence for v in self._v]
        )[:CDR3_OFFSET]

        for (start, size) in self.insertions:
            common_v = ''.join([
                common_v[:start],
                '.' * size,
                common_v[start:]
            ])

        common_v = common_v.replace('-', '')
        for gap in re.finditer('[-]+', self.germline):
            start, end = gap.span()
            common_v = ''.join([
                common_v[:start],
                '-' * (end - start),
                common_v[start:]
            ])

        self.germline = common_v
        self._v_end = len(self.germline)
        self._find_j(self._v_end)
        self.calculate_stats()

        offset = 0
        for i, region_size in enumerate(self.regions):
            for (start, size) in self.insertions:
                if offset <= start < offset + region_size:
                    self.regions[i] += size
            offset += region_size

        self.sequence = self.sequence.replace('.', '-')
        self.germline = self.germline.replace('.', '-')

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
            i = self.sequence[offset:].rfind(match)
            if i >= 0:
                return self._found_j(i + offset, j_gene, match, full_anchor)

            i = Seq(self.sequence).reverse_complement().rfind(match)
            if i >= 0:
                self.sequence = str(Seq(self.sequence).reverse_complement())
                if self.quality is not None:
                    self.quality.reverse()
                return self._found_j(i + offset, j_gene, match, full_anchor)
        raise AlignmentException('Could not find J anchor')

    def _found_j(self, i, j_gene, match, full_anchor):
        # If a match is found, record its location and gene
        self.j_anchor_pos = i
        self._j = [j_gene]

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

        self._j = self.j_germlines.get_ties(self.j_gene[0], match)
        self.j_length = len(j_full)

    def _find_v(self):
        for self._v_end in find_v_position(self.sequence):
            self._found_v()
            if self._v is not None:
                break

        if self._v is None:
            raise AlignmentException('Could not find suitable V anchor')

    def _found_v(self):
        '''Finds the V gene closest to that of the sequence'''
        if self._force_vs is not None:
            germlines = {vg: self.v_germlines[vg] for vg in self._force_vs}
        else:
            germlines = self.v_germlines

        aligned_v = VGene(self.sequence)
        v_score = None
        for v, germ in sorted(germlines.iteritems()):
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
                    self._germ_pos = germ.ungapped_anchor_pos
                    v_score = dist
                elif dist == v_score:
                    # Add the V-tie
                    self._v.append(v)

        if self._v is not None:
            # Determine the pad length
            self._pad_len = self._germ_pos - self._v_end
            # Mutation ratio is the distance divided by the length of overlap
            self.mutation_fraction = v_score / float(self.v_length)

    def align_to_germline(self, avg_len=None, avg_mut=None):
        if avg_len is not None and avg_mut is not None:
            self._v = self.v_germlines.get_ties(self.v_gene, avg_len, avg_mut)
        # Set the germline to the V gene up to the CDR3
        self.germline = get_common_seq(
            [self.v_germlines[v].sequence for v in self._v]
        )[:CDR3_OFFSET]
        self._pad_len = len(self.germline.replace('-', '')) - self._v_end
        # If we need to pad the sequence, do so, otherwise trim the sequence to
        # the germline length
        if self._pad_len >= 0:
            self.sequence = 'N' * self._pad_len + str(self.sequence)
            if self.quality is not None:
                self.quality = ([None] * self._pad_len) + self.quality
        else:
            self.removed_prefix = self.sequence[:-self._pad_len]
            self.sequence = str(self.sequence[-self._pad_len:])
            if self.quality is not None:
                self.removed_prefix_qual = self.quality[:-self._pad_len]
                self.quality = self.quality[-self._pad_len:]
        # Update the anchor positions after adding padding / trimming
        self.j_anchor_pos += self._pad_len
        self._v_end += self._pad_len

        # Add germline gaps to sequence before CDR3 and update anchor positions
        for i, c in enumerate(self.germline):
            if c == '-':
                self.sequence = self.sequence[:i] + '-' + self.sequence[i:]
                if self.quality is not None:
                    self.quality.insert(i, None)
                self.j_anchor_pos += 1
                self._v_end += 1

        self.calculate_stats()

    def calculate_stats(self):
        j_germ = get_common_seq(
            map(reversed, [self.j_germlines[j] for j in self.j_gene]))
        j_germ = ''.join(reversed(j_germ))
        # Calculate the length of the CDR3
        self._cdr3_len = (
            self.j_anchor_pos + self.j_germlines.anchor_len -
            self.j_germlines.upstream_of_cdr3 - self._v_end
        )

        if self._cdr3_len <= 0:
            self._v = None
            raise AlignmentException('CDR3 has non-positive length: {}'
                                     '.'.format(self._cdr3_len))

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
            # If the germline is longer than the sequence, there was probably a
            # deletion, so flag it as such
            self.sequence += 'N' * (len(self.germline) - len(self.sequence))
            if self.quality is not None:
                self.quality.extend([None] * (len(self.germline) -
                                     len(self.quality)))
            self._possible_indel = True

        # Get the pre-CDR3 germline and sequence stripped of gaps
        pre_cdr3_germ = self.germline[:self._v_end]
        pre_cdr3_seq = self.sequence[:self._v_end]
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
        post_s = self.sequence[self._v_end+len(self.cdr3):]

        # Calculate their match count
        self.post_cdr3_match = self.post_cdr3_length - dnautils.hamming(
            post_j, post_s)

        self.v_match = self.v_length - dnautils.hamming(
            self.germline[:self._v_end],
            self.sequence[:self._v_end]
        )

        self.j_match = self.j_length - dnautils.hamming(
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
        germ = self.germline[start:self._v_end]
        seq = self.sequence[start:self._v_end]

        for i in range(0, len(germ) - self.INDEL_WINDOW + 1):
            dist = dnautils.hamming(germ[i:i+self.INDEL_WINDOW],
                                    seq[i:i+self.INDEL_WINDOW])
            if dist >= self.INDEL_MISMATCH_THRESHOLD * self.INDEL_WINDOW:
                return True

        return False
