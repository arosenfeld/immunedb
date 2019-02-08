import re

from Bio.Seq import Seq

import dnautils
from immunedb.common.models import CDR3_OFFSET
import immunedb.util.lookups as lookups


class VDJSequence(object):
    def __init__(self, seq_id, sequence, quality=None, rev_comp=False,
                 copy_number=1):
        if quality and len(sequence) != len(quality):
            raise ValueError('Sequence and quality must be the same length')
        if not all([c in 'ATCGN-' for c in sequence]):
            raise ValueError('Invalid characters in sequence: {}'.format(
                sequence))

        self.seq_id = seq_id
        self.copy_number = copy_number
        self.orig_sequence = sequence[:]
        self.orig_quality = quality[:] if quality else None
        self.rev_comp = rev_comp
        self._sequence = sequence
        self._quality = quality
        self._removed_prefix_sequence = ''
        self._removed_prefix_quality = ''

    @property
    def sequence(self):
        return self._sequence

    @property
    def quality(self):
        return self._quality

    @property
    def removed_prefix_sequence(self):
        return self._removed_prefix_sequence

    @property
    def removed_prefix_quality(self):
        return self._removed_prefix_quality

    def reverse_complement(self):
        return VDJSequence(
            self.seq_id,
            str(Seq(self._sequence).reverse_complement()),
            self._quality[::-1] if self._quality else None,
            rev_comp=True,
            copy_number=self.copy_number
        )

    def pad(self, count):
        self._sequence = ('N' * count) + self._sequence
        if self._quality:
            self._quality = (' ' * count) + self._quality

    def pad_right(self, count):
        self._sequence += 'N' * count
        if self._quality:
            self._quality += ' ' * count

    def remove_prefix(self, count):
        self._removed_prefix_sequence = self._sequence[:count]
        self._sequence = self._sequence[count:]
        if self._quality:
            self._removed_prefix_quality = self._sequence[:count]
            self._quality = self._quality[count:]

    def trim(self, count):
        new_prefix = ''.join([
            c if c == '-' else 'N' for c in self._sequence[:count]
        ])
        self._sequence = new_prefix + self._sequence[count:]
        if self._quality:
            self._quality = (' ' * count) + self._quality[count:]

    def trim_right(self, count):
        self._sequence = self._sequence[:count]
        if self._quality:
            self._quality = self._quality[:count]

    def add_gap(self, pos, char='-'):
        self._sequence = self._sequence[:pos] + char + self._sequence[pos:]
        if self._quality:
            self._quality = self._quality[:pos] + ' ' + self._quality[pos:]

    def rfind(self, seq):
        return self._sequence.rfind(seq)

    def __getitem__(self, key):
        return self._sequence[key]

    def __setitem__(self, key, value):
        self._sequence[key] = value

    def __len__(self):
        return len(self._sequence)


class VDJAlignment(object):
    INDEL_WINDOW = 30
    INDEL_MISMATCH_THRESHOLD = .6

    def __init__(self, sequence):
        self.sequence = sequence
        self.germline = None
        self.j_gene = set()
        self.v_gene = set()

        self.locally_aligned = False
        self.seq_offset = 0
        self.v_length = 0
        self.j_length = 0
        self.v_mutation_fraction = 0
        self.cdr3_start = CDR3_OFFSET
        self.cdr3_num_nts = 0
        self.germline_cdr3 = None
        self.post_cdr3_length = 0
        self.insertions = set([])
        self.deletions = set([])

    @property
    def filled_germline(self):
        return ''.join((
            self.germline[:self.cdr3_start],
            self.germline_cdr3,
            self.germline[self.cdr3_start + self.cdr3_num_nts:]
        ))

    @property
    def seq_start(self):
        return max(0, self.seq_offset)

    @property
    def num_gaps(self):
        return self.sequence[self.seq_start:self.cdr3_start].count('-')

    @property
    def cdr3(self):
        return self.sequence[self.cdr3_start:self.cdr3_start +
                             self.cdr3_num_nts]

    @property
    def partial(self):
        return self.seq_start > 0

    @property
    def in_frame(self):
        return len(self.cdr3) % 3 == 0 and self.cdr3_start % 3 == 0

    @property
    def stop(self):
        return lookups.has_stop(self.sequence)

    @property
    def functional(self):
        return self.in_frame and not self.stop

    @property
    def v_match(self):
        start = self.seq_start
        end = start + self.v_length + self.num_gaps

        return self.v_length - dnautils.hamming(
            self.filled_germline[start:end],
            self.sequence[start:end]
        )

    @property
    def j_match(self):
        return self.j_length - dnautils.hamming(
            self.filled_germline[-self.j_length:],
            self.sequence[-self.j_length:]
        )

    @property
    def pre_cdr3_length(self):
        return self.cdr3_start - self.seq_start - self.num_gaps

    @property
    def pre_cdr3_match(self):
        start = self.seq_start + self.num_gaps
        end = self.cdr3_start

        return self.pre_cdr3_length - dnautils.hamming(
            self.germline[start:end],
            self.sequence[start:end]
        )

    @property
    def post_cdr3_match(self):
        return self.post_cdr3_length - dnautils.hamming(
            self.germline[-self.post_cdr3_length:],
            self.sequence[-self.post_cdr3_length:]
        )

    @property
    def has_possible_indel(self):
        # Start comparison on first full AA to the INDEL_WINDOW or CDR3,
        # whichever comes first
        start = re.search('[ATCG]', self.sequence.sequence).start()
        germ = self.germline[start:self.cdr3_start]
        seq = self.sequence[start:self.cdr3_start]

        for i in range(0, len(germ) - self.INDEL_WINDOW + 1):
            dist = dnautils.hamming(germ[i:i+self.INDEL_WINDOW],
                                    seq[i:i+self.INDEL_WINDOW])
            if dist >= self.INDEL_MISMATCH_THRESHOLD * self.INDEL_WINDOW:
                return True

        return False

    def trim_to(self, count):
        old_pad = self.seq_start - self.sequence[:self.seq_start].count('-')
        n_extension = re.match(
            '[N]*', self.sequence[count:]).span()[1]

        self.sequence.trim(count)
        self.seq_offset = re.match(r'[N-]*', self.sequence.sequence).span()[1]
        self.seq_offset -= n_extension

        new_pad = self.sequence[:self.seq_start].count('-')
        self.v_length = self.v_length - self.seq_start + old_pad + new_pad
