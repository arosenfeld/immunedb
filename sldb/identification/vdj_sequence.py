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
                self._find_v_position()
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
    def to_cdr3_length(self):
        return self._to_cdr3_length

    @property
    def to_cdr3_match(self):
        return self._to_cdr3_match

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

        for match, j_gene in anchors.all_j_anchors():
            i = seq.find(match)
            if i >= 0:
                if rev_comp:
                    self._j_anchor_pos = i
                else:
                    self._j_anchor_pos = len(self._seq) - i
                self._j = j_gene
                self._j_match = len(j_gene)
                self._j_length = len(j_gene)
                if rev_comp:
                    self._seq = self._seq.reverse_complement()
                return

    def _find_v_position(self):
        '''Tries to find the end of the V gene region'''
        # Try to find DxxxyzC
        dc = self._find_dc()
        if dc is not None:
            dc_seq = Seq(dc.group(0))
            # Make sure yz in ['YY', 'YC', 'YH']
            if str(dc_seq.translate()[-2:]) in anchors.dc_final_aas:
                self._v_anchor = dc_seq
                self._v_anchor_pos = dc.end() - 3
                return

        # If DxxyzC isn't found, try to find 'YYC', 'YCC', or 'YHC'
        yxc = self._find_yxc()
        if yxc is not None:
            self._v_anchor = Seq(yxc.group(0))
            self._v_anchor_pos = yxc.end() - 3

    def _find_v(self):
        '''Finds the V gene closest to that of the sequence'''
        self._v_score = None

        for v, germ in germlines.v.iteritems():
            # Strip the gaps
            v_seq = Seq(germ.replace('-', ''))
            # Get the last occurrence of 'TGT', the germline anchor
            germ_pos = v_seq.rfind(anchors.germline_anchor)
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
                self._to_cdr3_seq = s_seq[:self.v_anchor_pos]
                self._to_cdr3_germ = v_seq[:germ_pos - diff]
            elif dist == self._v_score:
                self._v.append(v)

        pad_len = germlines.v[self._v[0]].replace('-', '').rfind(
            anchors.germline_anchor) - self.v_anchor_pos

        # If we need to pad with a full sequence, there is a misalignment
        if self._full_v and pad_len > 0:
            self._v = None
            return

        self._to_cdr3_length = len(self._to_cdr3_seq)
        self._to_cdr3_match = self._to_cdr3_length - distance.hamming(
            str(self._to_cdr3_seq), str(self._to_cdr3_germ))

        self._germline = germlines.v[sorted(self._v)[0]][:self.CDR3_OFFSET]
        if pad_len >= 0:
            self._seq = 'N' * pad_len + str(self._seq)
        else:
            self._seq = str(self._seq[-1 * pad_len:])

        # Mutation ratio is the distance divided by the length of overlap
        self._mutation_frac = self._v_score / float(self._v_length)
        self._j_anchor_pos += pad_len
        self._v_anchor_pos += pad_len

        # Add germline gaps to sequence before CDR3 and update anchor positions
        for i, c in enumerate(self._germline):
            if c == '-':
                self._seq = self._seq[:i] + '-' + self._seq[i:]
                self._j_anchor_pos += 1
                self._v_anchor_pos += 1

        # Find the J anchor in the germline J gene
        j_anchor_in_germline = germlines.j[self.j_gene].rfind(
            str(anchors.j_anchors[self.j_gene]))
        # Calculate the length of the CDR3
        self._cdr3_len = self.j_anchor_pos - self.CDR3_OFFSET - \
            j_anchor_in_germline
        self._j_anchor_pos += self._cdr3_len
        # Fill germline CDR3 with gaps
        self._germline += '-' * self._cdr3_len
        self._germline += germlines.j[self.j_gene]
        self._seq = self._seq[:len(self._germline)]

    def _find_dc(self):
        '''Finds the first occurrence of the amino-acid sequence DxxxxxC'''
        for d in lookups.aa_to_all_nts('D'):
            for c in lookups.aa_to_all_nts('C'):
                find = re.search(d + '([ATCG]{15})' + c, str(self.sequence))
                if find is not None:
                    return find
        return None

    def _find_yxc(self):
        for aas in ['YYC', 'YCC', 'YHC']:
            for nts in lookups.aa_to_all_nts(aas):
                find = re.search(nts, str(self.sequence))
                if find is not None:
                    return find
        return None

class Best(object):
    def __init__(self, gene, score, match):
        self._genes = set([gene])
        self._score = score
        self._match = match

    @property
    def score(self):
        return self._score

    @property
    def genes(self):
        return self._genes

    @property
    def match(self):
        return self._match

    def add_gene(self, seq):
        self._genes.add(seq)
