import re
import distance

from Bio.Seq import Seq
import sldb.identification.anchors as anchors
import sldb.identification.germlines as germlines
import sldb.util.lookups as lookups


class VDJSequence(object):
    CDR3_OFFSET = 309
    MISMATCH_THRESHOLD = 3

    def __init__(self, id, seq, full_v):
        self._id = id
        self._seq = seq
        self._full_v = full_v
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

        # If there are Ns in the sequence, ignore it
        if 'N' not in self.sequence:
            self._find_j()
            if self._j is not None:
                self._v_anchor_pos = self._find_v_position(self.sequence)
                if self.v_anchor_pos is not None:
                    self._find_v()
                    if self.v_gene is not None:
                        self._calculate_stats()

    @property
    def id(self):
        return self._id

    @property
    def j_gene(self):
        return self._j

    @property
    def v_gene(self):
        return self._v

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
        self._compare(False)
        if self._j is None:
            self._compare(True)

    def _compare(self, rev_comp):
        '''Finds an exact J gene anchor match at the amino acid level'''
        if rev_comp:
            seq = self._seq.reverse_complement()
        else:
            seq = self._seq

        for match, full_anchor, j_gene in anchors.all_j_anchors():
            i = seq.rfind(match)
            if i >= 0:
                self._j_anchor_pos = i
                if rev_comp:
                    self._seq = self._seq.reverse_complement()
                self._j = j_gene

                j_full = germlines.j[self.j_gene]
                j_in_cdr3 = j_full[:len(j_full) - germlines.j_offset]
                cdr3_end = self._j_anchor_pos + len(match) - germlines.j_offset
                cdr3_segment = self.sequence[cdr3_end - len(j_in_cdr3):cdr3_end]
                j_full = j_full[self._find_streak_position(reversed(j_in_cdr3),
                                           reversed(cdr3_segment)):]
                streak = 0
                for i, (g, s) in enumerate(zip(reversed(j_in_cdr3),
                                               reversed(cdr3_segment))):
                    if g == s:
                        streak = 0
                    else:
                        streak += 1
                    if streak >= self.MISMATCH_THRESHOLD:
                        j_full = j_full[i:]
                        break

                self._j_start = self._j_anchor_pos + len(match) - len(j_full)

                if len(j_full) != len(
                        self.sequence[self._j_start:self._j_start+len(j_full)]):
                    self._j = None
                    self._j_anchor_pos = None
                    return

                dist = distance.hamming(
                    j_full,
                    self.sequence[self._j_start:self._j_start+len(j_full)])

                self._j_length = len(j_full)
                self._j_match = self._j_length - dist

                return

    def _find_v_position(self, sequence):
        '''Tries to find the end of the V gene region'''
        # Try to find DxxxyzC
        found = self._find_dc(sequence)
        if found is None:
            # If DxxyzC isn't found, try to find 'YYC', 'YCC', or 'YHC'
            found = self._find_yxc(sequence)
        return found

    def _find_v(self):
        '''Finds the V gene closest to that of the sequence'''
        self._v_score = None

        for v, germ in germlines.v.iteritems():
            # Strip the gaps
            v_seq = Seq(germ.replace('-', ''))
            # Find V anchor in germline
            germ_pos = self._find_v_position(v_seq)
            # Determine the gap between the two anchors
            diff = abs(germ_pos - self.v_anchor_pos)
            s_seq = self.sequence

            # Trim the sequence which has the maximal anchor position, and
            # determine the CDR3 start position without gaps
            if germ_pos > self.v_anchor_pos:
                v_seq = v_seq[diff:]
                cdr3_offset_in_v = self.CDR3_OFFSET - germ[diff:].count('-') - \
                    diff
            else:
                s_seq = s_seq[diff:]
                cdr3_offset_in_v = self.CDR3_OFFSET - germ.count('-')
            # Only compare to the start of J
            v_seq = v_seq[:self._j_start]
            s_seq = s_seq[:self._j_start]

            streak = 0
            v_cdr3 = v_seq[cdr3_offset_in_v:]
            s_cdr3 = s_seq[cdr3_offset_in_v:]
            max_index = cdr3_offset_in_v + (self._find_streak_position(
                v_cdr3, s_cdr3) - self.MISMATCH_THRESHOLD)
            v_seq = v_seq[:max_index]
            s_seq = s_seq[:max_index]

            # Determine the distance between the germline and sequence
            dist = distance.hamming(str(v_seq), str(s_seq))
            # Record this germline if it is has the lowest distance
            if self._v_score is None or dist < self._v_score:
                self._v = [v]
                self._v_score = dist
                self._v_length = len(s_seq)
                self._v_match = len(s_seq) - dist
                self._germ_pos = germ_pos
            elif dist == self._v_score:
                self._v.append(v)

        self._pad_len = self._germ_pos - self.v_anchor_pos

        # If we need to pad with a full sequence, there is a misalignment
        if self._full_v and self._pad_len > 0:
            self._v = None
            return

    def _calculate_stats(self):
        self._germline = germlines.v[sorted(self._v)[0]][:self.CDR3_OFFSET]
        if self._pad_len >= 0:
            self._seq = 'N' * self._pad_len + str(self._seq)
        else:
            self._seq = str(self._seq[-self._pad_len:])
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

        j_germ = germlines.j[self.j_gene]
        # Find the J anchor in the germline J gene
        j_anchor_in_germline = j_germ.rfind(
            str(anchors.j_anchors[self.j_gene]))
        # Calculate the length of the CDR3
        self._cdr3_len = (self.j_anchor_pos + \
            len(anchors.j_anchors[self.j_gene]) - germlines.j_offset) - \
            self.v_anchor_pos
        
        if self._cdr3_len <= 0:
            self._v = None
            return

        self._j_anchor_pos += self._cdr3_len
        # Fill germline CDR3 with gaps
        self._germline += '-' * self._cdr3_len
        self._germline += j_germ[-germlines.j_offset:]
        if len(self.sequence) > len(self.germline):
            self._seq = self._seq[:len(self._germline)]
            self.probable_deletion = False
        else:
            self._seq += 'N' * (len(self.germline) - len(self.sequence))
            self.probable_deletion = True

        self._num_gaps = self.sequence[:self.CDR3_OFFSET].count('-')

        pre_cdr3_germ = self.germline[:self.CDR3_OFFSET].replace('-', '')
        pre_cdr3_seq = self.sequence[:self.CDR3_OFFSET].replace('-', '')
        if self._pad_len > 0:
            pre_cdr3_germ = pre_cdr3_germ[self._pad_len:]
            pre_cdr3_seq = pre_cdr3_seq[self._pad_len:]

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

    def _find_dc(self, sequence):
        return self._find_with_frameshifts(sequence, 'D(.{3}((YY)|(YC)|(YH)))C')

    def _find_yxc(self, sequence):
        return self._find_with_frameshifts(sequence, 'Y([YHC])C')

    def _find_with_frameshifts(self, sequence, regex):
        for shift in range(0, 3):
            seq = sequence[shift:]
            seq = seq[:len(seq) - len(seq) % 3]
            aas = str(seq.translate())
            res = re.search(regex, aas)
            if res is not None:
                return (res.end() - 1) * 3 + shift
        return None

    def _find_streak_position(self, s1, s2):
        streak = 0
        for i, (c1, c2) in enumerate(zip(s1, s2)):
            streak = streak + 1 if c1 != c2 else 0
            if streak >= self.MISMATCH_THRESHOLD:
                break
        return i
