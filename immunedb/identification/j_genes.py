import dnautils
from Bio import SeqIO
from immunedb.identification import GeneTies


class JGermlines(GeneTies):
    ALLOWED_PREFIX = set([
        'IGHJ',
        'TRAJ',
        'TRBJ',
        'TRDJ',
        'TRGJ'
    ])
    def __init__(self, path_to_germlines, upstream_of_cdr3, anchor_len,
                 min_anchor_len, ties_prob_threshold=.01):
        self._upstream_of_cdr3 = upstream_of_cdr3
        self._anchor_len = anchor_len
        self._min_anchor_len = min_anchor_len
        self._min_length = None

        with open(path_to_germlines) as fh:
            for record in SeqIO.parse(fh, 'fasta'):
                self.prefix = record.id[:4]
                assert self.prefix in JGermlines.ALLOWED_PREFIX
                if all(map(lambda c: c in 'ATCGN-', record.seq)):
                    self[record.id] = str(record.seq).upper()
                    if (self._min_length is None or
                            len(self[record.id]) < self._min_length):
                        self._min_length = len(self[record.id])

        self._anchors = {name: seq[-anchor_len:] for name, seq in
                         self.iteritems()}
        super(JGermlines, self).__init__(
            {k: v for k, v in self.iteritems()},
            ties_prob_threshold=ties_prob_threshold
        )

    @property
    def upstream_of_cdr3(self):
        return self._upstream_of_cdr3

    @property
    def anchor_len(self):
        return self._anchor_len

    @property
    def full_anchors(self):
        return self._anchors

    def get_j_in_cdr3(self, gene):
        return self[gene][:-self._upstream_of_cdr3]

    def get_all_anchors(self, allowed_genes=None):
        if allowed_genes is None:
            allowed_genes = self
        else:
            allowed_genes = {k: v for k, v in self.iteritems() if k in
                             allowed_genes}
        max_len = max(map(len, allowed_genes.values()))
        for trim_len in range(0, max_len, 3):
            for j, seq in allowed_genes.iteritems():
                trimmed_seq = seq[-self.anchor_len:-trim_len]
                if len(trimmed_seq) >= self._min_anchor_len:
                    yield trimmed_seq, j

    def get_single_tie(self, gene, length, mutation):
        seq = self[gene][-self.anchor_len:]
        tied = self.all_alleles(set([gene]))
        for j, other_seq in sorted(self.iteritems()):
            other_seq = other_seq[-self.anchor_len:][:len(seq)]
            if other_seq == seq:
                tied.add(j)
            elif dnautils.hamming(other_seq, seq) == 0:
                tied.add(j)
        return tied
