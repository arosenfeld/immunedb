import re
import distance

from Bio.Seq import Seq
import sldb.identification.anchors as anchors
import sldb.identification.germlines as germlines
import sldb.util.lookups as lookups

def find_v_position(sequence, reverse):
    '''Tries to find the end of the V gene region'''
    loc = sequence.rfind('TATTACTGT')
    if loc >= 0:
        return loc + 6
    # Try to find DxxxyzC
    found = find_dc(sequence, reverse)
    if found is None:
        # If DxxyzC isn't found, try to find 'YYC', 'YCC', or 'YHC'
        found = find_yxc(sequence, reverse)
    return found

def find_dc(sequence, reverse):
    return find_with_frameshifts(sequence, 'D(.{3}((YY)|(YC)|(YH)))C', reverse)

def find_yxc(sequence, reverse):
    return find_with_frameshifts(sequence, 'Y([YHC])C', reverse)

def find_with_frameshifts(sequence, regex, reverse):
    r = range(2, -1, -1)
    if reverse:
        r = reversed(r)

    for shift in r:
        seq = sequence[shift:]
        seq = seq[:len(seq) - len(seq) % 3]
        aas = str(seq.translate())
        res = re.search(regex, aas)
        if res is not None:
            return (res.end() - 1) * 3 + shift
    return None


class VDJSequence(object):
    CDR3_OFFSET = 309
    MISMATCH_THRESHOLD = 3
    MIN_J_ANCHOR_LEN = 12

    def __init__(self, id, seq, full_v, v_germlines):
        self._id = id
        self._seq = seq
        self._full_v = full_v
        self.v_germlines = v_germlines

        self._seq_filled = None

        self._j = None
        self._j_anchor_pos = None
        self._j_match = None

        self._v = None
        self._v_anchor_pos = None
        self._v_match = None

        self._mutation_frac = None
        self._germline = None
        self._cdr3_len = 0

        self.probable_deletion = False

        # If there are invalid characters in the sequence, ignore it 
        stripped = filter(lambda s: s in 'ATCGN', self.sequence)
        if len(stripped) == len(self.sequence):
            self._find_j()
            if self._j is not None:
                self._get_v(reverse=False)


    @property
    def id(self):
        return self._id

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
        return self.sequence[self.CDR3_OFFSET:
                             self.CDR3_OFFSET + self._cdr3_len]

    @property
    def sequence(self):
        return self._seq

    @property
    def sequence_filled(self):
        if self._seq_filled is None:
            self._seq_filled = ''
            for i, c in enumerate(self.sequence):
                if c.upper() == 'N':
                    self._seq_filled += self.germline[i].upper()
                else:
                    self._seq_filled += c
        return self._seq_filled

    @property
    def germline(self):
        return self._germline

    @property
    def mutation_fraction(self):
        return self._mutation_frac

    @property
    def num_gaps(self):
        return self._num_gaps

    @property
    def pad_length(self):
        return self._pad_len if self._pad_len >=0 else 0

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

    def _find_j(self):
        '''Finds the location and type of J gene'''
        # Iterate over every possible J anchor.  For each germline, try its full
        # sequence, then exclude the final 3 characters at a time until
        # there are only MIN_J_ANCHOR_LEN nucleotides remaining.
        #
        # For example, the order for one germline:
        # TGGTCACCGTCTCCTCAG
        # TGGTCACCGTCTCCT
        # TGGTCACCGTCT
        for match, full_anchor, j_gene in anchors.all_j_anchors(
                self.MIN_J_ANCHOR_LEN):
            i = self._seq.rfind(match)
            if i >= 0:
                return self._found_j(i, j_gene, match, full_anchor)

            i = self._seq.reverse_complement().rfind(match)
            if i >= 0:
                # If found in the reverse complement, flip and translate the
                # actual sequence for the rest of the analysis
                self._seq = self._seq.reverse_complement()
                return self._found_j(i, j_gene, match, full_anchor)

    def _found_j(self, i, j_gene, match, full_anchor):
        # If a match is found, record its location and gene
        self._j_anchor_pos = i
        self._j = [j_gene]

        # Get the full germline J gene
        j_full = germlines.j[self.j_gene[0]]
        # Get the portion of J in the CDR3
        j_in_cdr3 = j_full[:len(j_full) - germlines.j_offset]
        cdr3_end = (self._j_anchor_pos) - germlines.j_offset +\
            len(match)
        cdr3_segment = self.sequence[cdr3_end - len(j_in_cdr3):cdr3_end]
        if len(j_in_cdr3) == 0 or len(cdr3_segment) == 0:
            self._j = None
            return
        # Get the extent of the J in the CDR3
        streak = self._find_streak_position(reversed(j_in_cdr3),
                                   reversed(cdr3_segment))
        # Trim the J gene based on the extent in the CDR3
        if streak is not None:
            j_full = j_full[streak:]

        # Find where the full J starts
        self._j_start = self._j_anchor_pos + len(match) - len(j_full)

        # If the trimmed germline J extends past the end of the
        # sequence, there is a misalignment
        if len(j_full) != len(
                self.sequence[self._j_start:self._j_start+len(j_full)]):
            self._j = None
            self._j_anchor_pos = None
            return

        # Get the full-J distance
        dist = distance.hamming(
            j_full,
            self.sequence[self._j_start:self._j_start+len(j_full)])

        self._j = anchors.get_j_ties(self.j_gene[0], match)
        self._j_length = len(j_full)
        self._j_match = self._j_length - dist

    def _get_v(self, reverse):
        self._v_anchor_pos = find_v_position(self.sequence, reverse)
        if self.v_anchor_pos is not None:
            self._find_v()
            if self.v_gene is not None:
                self._calculate_stats()
            elif not reverse:
                self._get_v(True)

    def _find_v(self):
        '''Finds the V gene closest to that of the sequence'''
        self._v_score = None
        self._v = None

        for v, (germ, germ_pos) in sorted(self.v_germlines.iteritems()):
            dist, s_seq = self._compare_to_germline(germ, germ_pos)
            # Record this germline if it is has the lowest distance
            if dist is not None:
                if self._v_score is None or dist < self._v_score:
                    self._v = [v]
                    self._v_score = dist
                    self._v_length = len(s_seq)
                    self._v_match = len(s_seq) - dist
                    self._germ_pos = germ_pos
                elif dist == self._v_score:
                    # Add the V-tie
                    self._v.append(v)

        if self._v is None:
            return
        # Determine the pad length
        self._pad_len = self._germ_pos - self.v_anchor_pos

        # If we need to pad with a full sequence, there is a misalignment
        if self._full_v and self._pad_len > 0:
            self._v = None
            return
        self._v = list(set(map(lambda v: v.split('*')[0], self._v)))

    def _compare_to_germline(self, germ, germ_pos):
        # Strip the gaps
        v_seq = Seq(germ.replace('-', ''))
        # Find V anchor in germline
        if germ_pos is None:
            return None, None
        # Determine the gap between the two anchors
        diff = abs(germ_pos - self.v_anchor_pos)
        s_seq = self.sequence

        # Trim the sequence which has the maximal anchor position, and
        # determine the CDR3 start position without gaps
        if germ_pos > self.v_anchor_pos:
            v_seq = v_seq[diff:]
            cdr3_offset_in_v = germ_pos - diff
        else:
            s_seq = s_seq[diff:]
            cdr3_offset_in_v = germ_pos

        # Only compare to the start of J
        v_seq = v_seq[:self._j_start]
        s_seq = s_seq[:self._j_start]

        # Determine the CDR3 in the germline and sequence
        v_cdr3 = v_seq[cdr3_offset_in_v:]
        s_cdr3 = s_seq[cdr3_offset_in_v:]
        v_cdr3 = v_cdr3[:min(len(v_cdr3), len(s_cdr3))]
        s_cdr3 = s_cdr3[:min(len(v_cdr3), len(s_cdr3))]
        if len(v_cdr3) == 0 or len(s_cdr3) == 0:
            return None, None

        # Find the extent of the sequence's V into the CDR3
        streak = self._find_streak_position(v_cdr3, s_cdr3)
        if streak is not None:
            # If there is a streak of mismatches, cut after the streak
            max_index = cdr3_offset_in_v + (streak - self.MISMATCH_THRESHOLD)
        else:
            # Unlikely: the CDR3 in the sequence exactly matches the
            # germline.  Use the smaller sequence length (full match)
            max_index = cdr3_offset_in_v + min(len(v_cdr3), len(s_cdr3))
        # Compare to the end of V
        v_seq = v_seq[:max_index]
        s_seq = s_seq[:max_index]

        # Determine the distance between the germline and sequence
        dist = distance.hamming(str(v_seq), str(s_seq))

        return dist, s_seq

    def _calculate_stats(self):
        # Set the germline to the V gene up to the CDR3
        self._germline = self.v_germlines[sorted(self._v)[0] +
            '*01'][0][:self.CDR3_OFFSET]
        # If we need to pad the sequence, do so, otherwise trim the sequence to
        # the germline length
        if self._pad_len >= 0:
            self._seq = 'N' * self._pad_len + str(self._seq)
        else:
            self._seq = str(self._seq[-self._pad_len:])
        # Update the anchor positions after adding padding / trimming
        self._j_anchor_pos += self._pad_len
        self._v_anchor_pos += self._pad_len

        # Mutation ratio is the distance divided by the length of overlap
        self._mutation_frac = self._v_score / float(self._v_length)

        # Add germline gaps to sequence before CDR3 and update anchor positions
        for i, c in enumerate(self._germline):
            if c == '-':
                self._seq = self._seq[:i] + '-' + self._seq[i:]
                self._j_anchor_pos += 1
                self._v_anchor_pos += 1

        j_germ = germlines.j[self.j_gene[0]]
        # Find the J anchor in the germline J gene
        j_anchor_in_germline = j_germ.rfind(
            str(anchors.j_anchors[self.j_gene[0]]))
        # Calculate the length of the CDR3
        self._cdr3_len = (self.j_anchor_pos + \
            len(anchors.j_anchors[self.j_gene[0]]) - germlines.j_offset) - \
            self.v_anchor_pos
        
        if self._cdr3_len <= 0:
            self._v = None
            return

        self._j_anchor_pos += self._cdr3_len
        # Fill germline CDR3 with gaps
        self._germline += '-' * self._cdr3_len
        self._germline += j_germ[-germlines.j_offset:]
        # If the sequence is longer than the germline, trim it
        if len(self.sequence) > len(self.germline):
            self._seq = self._seq[:len(self._germline)]
        elif len(self.sequence) < len(self.germline):
            # If the germline is longer than the sequence, there was probably a
            # deletion, so flag it as such
            self._seq += 'N' * (len(self.germline) - len(self.sequence))
            self.probable_deletion = True

        # Get the number of gaps
        self._num_gaps = self.sequence[:self.CDR3_OFFSET].count('-')

        # Get the pre-CDR3 germline and sequence stripped of gaps
        pre_cdr3_germ = self.germline[:self.CDR3_OFFSET].replace('-', '')
        pre_cdr3_seq = self.sequence[:self.CDR3_OFFSET].replace('-', '')
        # If there is padding, get rid of it in the sequence and align the
        # germline
        if self._pad_len > 0:
            pre_cdr3_germ = pre_cdr3_germ[self._pad_len:]
            pre_cdr3_seq = pre_cdr3_seq[self._pad_len:]

        # Calculate the pre-CDR3 length and distance
        self._pre_cdr3_length = len(pre_cdr3_seq)
        self._pre_cdr3_match = self._pre_cdr3_length - distance.hamming(
            str(pre_cdr3_seq), str(pre_cdr3_germ))

        # Get the length of J after the CDR3
        self._post_cdr3_length = germlines.j_offset
        # Get the sequence and germline sequences after CDR3
        post_j = j_germ[-self.post_cdr3_length:]
        post_s = self.sequence[self.CDR3_OFFSET+len(self.cdr3):]

        # Calculate their match count
        self._post_cdr3_match = self.post_cdr3_length - distance.hamming(
            post_j, post_s)

    def _find_streak_position(self, s1, s2):
        ''' Finds the first streak of MISMATCH_THRESHOLD characters where s1
        does not equal s2.

        For example if MISMATCH_THRESHOLD is 3:
            ATCGATCGATCGATCG
            ATCGATCGATCTTACG
                         ^--- Returned index
        '''
        streak = 0
        for i, (c1, c2) in enumerate(zip(s1, s2)):
            streak = streak + 1 if c1 != c2 else 0
            if streak >= self.MISMATCH_THRESHOLD:
                return i
        return None
