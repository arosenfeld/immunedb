import re

from Bio.Seq import Seq
import sldb.identification.anchors as anchors
import sldb.identification.v_genes as v_genes 
import sldb.util.lookups as lookups


def _hamming_distance(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


class VDJSequence(object):
    def __init__(self, seq):
        self._seq = seq
        self._j = None
        self._j_anchor_pos = None
        self._v = None
        self._v_anchor_pos = None
        self._is_rev = False
        self._mutation_frac = None

        self._find_j()
        if self._j is not None:
            self._find_v_position()
            if self.v_anchor_pos is not None:
                self._find_v()

    @property
    def j_gene(self):
        return self._j

    @property
    def v_gene(self):
        return self._v

    @property
    def j_anchor_pos(self):
        return self._j_anchor_pos

    @property
    def v_anchor_pos(self):
        return self._v_anchor_pos

    @property
    def sequence(self):
        if self._is_rev:
            return self._seq.reverse_complement()
        return self._seq

    @property
    def mutation_fraction(self):
        return self._mutation_frac

    def _find_j(self):
        self._compare(False)
        if self._j is None:
            self._compare(True)

    def _compare(self, rev_comp):
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
                self._is_rev = rev_comp
                return

    def _find_v_position(self):
        # Try to find DXXXYZC
        dc = self._find_dc()
        if dc is not None:
            dc_seq = Seq(dc.group(0))
            # Make sure YZ in ['YY', 'YC', 'YH']
            if str(dc_seq.translate()[-2:]) in anchors.dc_final_aas:
                self._v_anchor = dc_seq
                self._v_anchor_pos = dc.end() - 2
                return

        # If DXXXYZC isn't found, try to find 'YYC', 'YCC', or 'YHC'
        yxc = self._find_yxc()
        if yxc is not None:
            self._v_anchor = Seq(yxc.group(0))
            self._v_anchor_pos = yxc.end() - 2
            return

    def _find_v(self):
        v_best = []
        v_score = None
        for v, germ in v_genes.germlines.iteritems():
            germ = Seq(germ)
            germ_pos = germ.rfind(anchors.germline_anchor)
            v_seq = germ[germ_pos - self.v_anchor_pos + 1:germ_pos + 1]
            s_seq = self.sequence[:self.v_anchor_pos]
            dist = _hamming_distance(str(v_seq), str(s_seq))
            if v_score is None or dist < v_score:
                v_best = [v]
                v_score = dist
            elif dist == v_score:
                v_best.append(v)

        self._v = v_best
        #TODO: is this padding correct for multiple v_bests
        pad = v_genes.germlines[v_best[0]].rfind(
            anchors.germline_anchor)
        self._seq = 'N' * pad + self._seq
        self._mutation_frac = v_score / float(self.v_anchor_pos)

    def _find_dc(self):
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
