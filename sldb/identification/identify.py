import sys
import re
from Levenshtein import distance as lev_distance

from Bio.Seq import Seq
from Bio import SeqIO

import sldb.identification.anchors as anchors
import sldb.util.lookups as lookups

def _hamming_distance(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


class VDJSequence(object):
    def __init__(self, seq):
        self._seq = seq
        self._j = None
        self._j_pos = None
        self._v = None
        self._v_anchor_pos = None
        self._is_rev = False

        self._find_j()
        if self._j is not None:
            self._find_v()

    @property
    def j_gene(self):
        return self._j

    @property
    def v_gene(self):
        return self._v

    @property
    def j_pos(self):
        return self._j_pos

    @property
    def v_anchor_pos(self):
        return self._v_anchor_pos

    @property
    def sequence(self):
        if self._is_rev:
            return self._seq.reverse_complement()
        return self._seq

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
                    self._j_pos = i
                else:
                    self._j_pos = len(self._seq) - i
                self._j = j_gene
                self._is_rev = rev_comp
                return

    def _find_v(self):
        dc = self._find_dc()
        if dc is not None:
            dc_seq = Seq(dc.group(1))
            if str(dc_seq.translate()[-2:]) in anchors.dc_final_aas:
                self._v_anchor = dc_seq
                self._v_anchor_pos = dc.end() - 1

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


if __name__ == '__main__':
    fn = 'fasta/D149-Ileum1-A703-A503-SEQ-14081-Wenzhao_S27_L001_R1_001.trim.fasta'
    with open(fn) as fh:
        for record in SeqIO.parse(fn, 'fasta'):
            vdj = VDJSequence(record.seq)
            if vdj.v_anchor_pos is not None and vdj.j_gene is not None:
                print 'found one'
                with open('v_gene.fasta') as fasta_h:
                    v_best = []
                    v_score = None
                    for v_record in SeqIO.parse(fasta_h, 'fasta'):
                        offset = vdj.v_anchor_pos
                        dist = _hamming_distance(str(v_record.seq[-offset:]),
                                                 str(vdj.sequence[:offset]))
                        #print dist, v_record.description
                        print v_record.seq[-offset:], v_record.description, dist
                        if v_score is None:
                            v_best = [v_record.description]
                            v_score = dist
                        else:
                            if dist == v_score:
                                v_best.append(v_record.description)
                            elif dist < v_score:
                                v_score = dist
                                v_best = [v_record.description]

                    print vdj.sequence[:offset]
                    print record.description
                    print vdj.sequence
                    print ' '*vdj.v_anchor_pos + vdj._v_anchor
                    print v_best
                    print v_score
                    break
