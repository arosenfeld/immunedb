import re
import distance

from Bio.Seq import Seq
import sldb.identification.anchors as anchors
import sldb.identification.germlines as germlines
import sldb.util.lookups as lookups


class VDJSequence(object):
    CDR3_OFFSET = 309

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
                j_start = self._j_anchor_pos - (len(j_full) - len(match))

                if len(j_full) != len(
                        self.sequence[j_start:j_start+len(j_full)]):
                    self._j = None
                    self._j_anchor_pos = None
                    return

                dist = distance.hamming(
                    j_full,
                    self.sequence[j_start:j_start+len(j_full)])

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
            germ_pos = self._find_v_position(v_seq)
            # Align the sequences for comparison
            diff = abs(germ_pos - self.v_anchor_pos)
            s_seq = self.sequence
            if germ_pos > self.v_anchor_pos:
                v_seq = v_seq[diff:]
            else:
                s_seq = s_seq[diff:]
            s_seq = s_seq[:len(v_seq)]

            # Determine the distance between the germline and sequence
            dist = distance.hamming(str(v_seq), str(s_seq))
            # Record this germline if it is has the lowest distance
            if self._v_score is None or dist < self._v_score:
                self._v = [v]
                self._v_score = dist
                self._v_length = len(s_seq)
                self._v_match = len(s_seq) - dist
                self._germ_pos = germ_pos
                if germ_pos > self.v_anchor_pos:
                    self._pre_cdr3_seq = s_seq[:self.v_anchor_pos]
                    self._pre_cdr3_germ = v_seq[:germ_pos - diff]
                else:
                    self._pre_cdr3_seq = s_seq[:self.v_anchor_pos - diff]
                    self._pre_cdr3_germ = v_seq[:germ_pos]

            elif dist == self._v_score:
                self._v.append(v)

        pad_len = self._germ_pos - self.v_anchor_pos

        # If we need to pad with a full sequence, there is a misalignment
        if self._full_v and pad_len > 0:
            self._v = None
            return

        self._pre_cdr3_length = len(self._pre_cdr3_seq)
        self._pre_cdr3_match = self._pre_cdr3_length - distance.hamming(
            str(self._pre_cdr3_seq), str(self._pre_cdr3_germ))

        self._germline = germlines.v[sorted(self._v)[0]][:self.CDR3_OFFSET]
        if pad_len >= 0:
            self._seq = 'N' * pad_len + str(self._seq)
        else:
            self._seq = str(self._seq[-pad_len:])
        self._j_anchor_pos += pad_len
        self._v_anchor_pos += pad_len

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
