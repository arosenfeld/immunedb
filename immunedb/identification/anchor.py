import re

import dnautils

from immunedb.common.models import CDR3_OFFSET
from immunedb.identification import AlignmentException, get_common_seq
from immunedb.identification.genes import VGene, find_v_position
from immunedb.identification.vdj_sequence import VDJAlignment


def sliding_window_match(sequence, match):
    r = re.search(match.replace('N', '.'), sequence)
    return r.start() if r else -1


class AnchorAligner(object):
    MISMATCH_THRESHOLD = 3

    def __init__(self, v_germlines, j_germlines):
        self.v_germlines = v_germlines
        self.j_germlines = j_germlines

    def get_alignment(self, vdj_sequence, limit_vs=None, limit_js=None):
        alignment = VDJAlignment(vdj_sequence)
        self.find_j(alignment, limit_js)
        self.find_v(alignment, limit_vs)
        return alignment

    def _find_index(self, sequence, germline):
        best_pos, best_hamming = None, None
        for pos in range(len(sequence) - len(germline)):
            hamming = dnautils.hamming(sequence[pos:pos + len(germline)],
                                       germline) / len(germline)
            if best_hamming is None or hamming < best_hamming:
                best_pos = pos
                best_hamming = hamming
                is_rc = False

        rc = sequence.reverse_complement()
        for pos in range(len(rc) - len(germline) + 1):
            hamming = dnautils.hamming(rc[pos:pos + len(germline)],
                                       germline) / len(germline)
            if best_hamming is None or hamming < best_hamming:
                best_pos = pos
                best_hamming = hamming
                is_rc = True

        best_pos += len(germline) - self.j_germlines.anchor_len
        return best_pos, best_hamming, is_rc

    def find_j(self, alignment, limit_js):
        # Iterate over every possible J anchor.  For each germline, try its
        # full sequence, then exclude the final 3 characters at a time until
        # there are only MIN_J_ANCHOR_LEN nucleotides remaining.
        #
        # For example, the order for one germline:
        # TGGTCACCGTCTCCTCAG
        # TGGTCACCGTCTCCT
        # TGGTCACCGTCT

        rc = alignment.sequence.reverse_complement()
        for match, j_gene in self.j_germlines.get_all_anchors(limit_js):
            i = alignment.sequence.rfind(match)
            if i >= 0:
                return self.process_j(alignment, i, len(match), limit_js)

            i = rc.rfind(match)
            if i >= 0:
                alignment.sequence = rc
                return self.process_j(alignment, i, len(match), limit_js)

            i = sliding_window_match(alignment.sequence.sequence, match)
            if i >= 0:
                return self.process_j(alignment, i, len(match), limit_js)

            i = sliding_window_match(rc.sequence, match)
            if i >= 0:
                alignment.sequence = rc
                return self.process_j(alignment, i, len(match), limit_js)

        # Last chance, find any matching position
        total_best_hamming = None
        total_best_rc, total_best_pos = None, None
        for j_gene, germ_seq in self.j_germlines.items():
            best_pos, best_hamming, is_rc = self._find_index(
                alignment.sequence, germ_seq)
            if (total_best_hamming is None or
                    best_hamming < total_best_hamming):
                total_best_hamming = best_hamming
                total_best_pos = best_pos
                total_best_rc = is_rc

        if total_best_rc:
            alignment.sequence = rc

        return self.process_j(alignment, total_best_pos,
                              self.j_germlines.anchor_len, limit_js)

    def process_j(self, alignment, i, match_len, limit_js):
        # If a match is found, record its location and gene
        alignment.j_anchor_pos = i
        end_of_j = min(
            alignment.j_anchor_pos + self.j_germlines.anchor_len,
            len(alignment.sequence)
        )
        if limit_js:
            j_germs = {
                k: v for k, v in self.j_germlines.items()
                if k.name in limit_js
            }
        else:
            j_germs = self.j_germlines

        best_dist = None
        for j_gene, j_seq in j_germs.items():
            seq_j = alignment.sequence[end_of_j - len(j_seq):end_of_j]
            if seq_j:
                dist = dnautils.hamming(seq_j, j_seq[:len(seq_j)]) / len(seq_j)
                if best_dist is None or dist < best_dist:
                    best_dist = dist
                    alignment.j_gene = set([j_gene])
                elif dist == best_dist:
                    alignment.j_gene.add(j_gene)

        if len(alignment.j_gene) == 0:
            raise AlignmentException('Could not find suitable J anchor')

        # Get the full germline J gene
        ex_j = sorted(alignment.j_gene)[0]
        j_full = self.j_germlines[ex_j]

        # Get the portion of the germline J in the CDR3
        germline_in_cdr3 = self.j_germlines.get_j_in_cdr3(ex_j)
        cdr3_end_pos = (
            alignment.j_anchor_pos + self.j_germlines.anchor_len -
            self.j_germlines.upstream_of_cdr3
        )
        sequence_in_cdr3 = alignment.sequence[
            cdr3_end_pos - len(germline_in_cdr3):
            cdr3_end_pos
        ]
        if len(germline_in_cdr3) == 0 or len(sequence_in_cdr3) == 0:
            alignment.j_gene = set()
            raise AlignmentException('Could not find sequence or germline in '
                                     'CDR3')

        # Get the extent of the J in the CDR3
        streak = dnautils.find_streak_position(
            germline_in_cdr3[::-1],
            sequence_in_cdr3[::-1],
            self.MISMATCH_THRESHOLD)

        # Trim the J gene based on the extent in the CDR3
        if streak is not None:
            j_full = j_full[len(germline_in_cdr3) - streak:]
            alignment.germline_cdr3 = germline_in_cdr3[-streak:]
        else:
            alignment.germline_cdr3 = germline_in_cdr3

        # Find where the full J starts
        j_start = alignment.j_anchor_pos + match_len - len(j_full)

        # If the trimmed germline J extends past the end of the
        # sequence, there is a misalignment
        if len(j_full) != len(alignment.sequence[j_start:j_start+len(j_full)]):
            alignment.j_gene = set()
            raise AlignmentException('Germline extended past end of J')

        alignment.j_length = len(j_full)
        alignment.post_cdr3_length = self.j_germlines.upstream_of_cdr3

    def find_v(self, alignment, limit_vs):
        for anchor_pos in find_v_position(alignment.sequence.sequence):
            self.process_v(alignment, anchor_pos, limit_vs)
            if len(alignment.v_gene) > 0:
                break

        if len(alignment.v_gene) == 0:
            raise AlignmentException('Could not find suitable V anchor')

    def process_v(self, alignment, anchor_pos, limit_vs):
        aligned_v = VGene(alignment.sequence.sequence)
        v_score = None
        for v, germ in self.v_germlines.alignments.items():
            if limit_vs is not None and v.name not in limit_vs:
                continue
            try:
                dist, total_length = germ.compare(aligned_v,
                                                  alignment.j_anchor_pos,
                                                  self.MISMATCH_THRESHOLD)
            except Exception:
                continue
            # Record this germline if it is has the lowest distance
            if dist is not None:
                if v_score is None or dist < v_score:
                    alignment.v_gene = set([v])
                    alignment.v_length = total_length
                    germ_pos = germ.ungapped_anchor_pos
                    v_score = dist
                elif dist == v_score:
                    # Add the V-tie
                    alignment.v_gene.add(v)

        if len(alignment.v_gene) > 0:
            # Determine the pad length
            alignment.seq_offset = germ_pos - anchor_pos
            # Mutation ratio is the distance divided by the length of overlap
            alignment.v_mutation_fraction = v_score / alignment.v_length

    def align_to_germline(self, alignment, avg_len=None, avg_mut=None):
        if avg_len is not None and avg_mut is not None:
            alignment.v_gene = self.v_germlines.get_ties(
                alignment.v_gene, avg_len, avg_mut)
            alignment.j_gene = self.j_germlines.get_ties(
                alignment.j_gene, avg_len, avg_mut)
        # Set the germline to the V gene up to the CDR3
        germ = get_common_seq(
            [self.v_germlines[v] for v in alignment.v_gene], cutoff=False
        )
        alignment.germline = germ[:CDR3_OFFSET]
        # If we need to pad the sequence, do so, otherwise trim the sequence to
        # the germline length
        if alignment.seq_offset >= 0:
            alignment.sequence.pad(alignment.seq_offset)
        else:
            alignment.sequence.remove_prefix(-alignment.seq_offset)
        alignment.j_anchor_pos += alignment.seq_offset

        # Add germline gaps to sequence before CDR3 and update anchor positions
        for i, c in enumerate(alignment.germline):
            if c == '-':
                alignment.sequence.add_gap(i)
                alignment.j_anchor_pos += 1
                if i < alignment.seq_start:
                    alignment.seq_offset += 1

        j_germ = get_common_seq(
            [self.j_germlines[j] for j in alignment.j_gene], right=True
        )
        # Calculate the length of the CDR3
        alignment.cdr3_num_nts = (
            alignment.j_anchor_pos + self.j_germlines.anchor_len -
            self.j_germlines.upstream_of_cdr3 - alignment.cdr3_start
        )

        v_end = alignment.seq_start + alignment.num_gaps + alignment.v_length
        v_germ = germ[CDR3_OFFSET:v_end]
        alignment.germline_cdr3 = ''.join((
            v_germ,
            '-' * (alignment.cdr3_num_nts -
                   len(v_germ) -
                   len(alignment.germline_cdr3)),
            alignment.germline_cdr3
        ))

        if alignment.cdr3_num_nts < 3:
            raise AlignmentException('CDR3 has no AAs')

        alignment.j_anchor_pos += alignment.cdr3_num_nts
        # Fill germline CDR3 with gaps
        alignment.germline += '-' * alignment.cdr3_num_nts
        alignment.germline += j_germ[-self.j_germlines.upstream_of_cdr3:]
        # If the sequence is longer than the germline, trim it
        if len(alignment.sequence) > len(alignment.germline):
            alignment.sequence.trim_right(len(alignment.germline))
        elif len(alignment.sequence) < len(alignment.germline):
            alignment.sequence.pad_right(len(alignment.germline) -
                                         len(alignment.sequence))
