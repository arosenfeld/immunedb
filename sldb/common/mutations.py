import sldb.util.lookups as lookups


class MutationType(object):
    """An enum-like class for different mutation types."""
    MUT_UNK = ('?', 'unknown')
    MUT_SYN = ('S', 'synonymous')
    MUT_CONS = ('C', 'conservative')
    MUT_UNCONS = ('U', 'nonconservative')
    MUT_NONE = (' ', 'none')

    @classmethod
    def get_symbol(cls, mtype):
        """Gets the single-character symbol for a mutation"""
        return mtype[0]

    @classmethod
    def get_readable(cls, mtype):
        """Gets the readable text for a mutation"""
        return mtype[1]

    @classmethod
    def get_types(cls):
        """Enumerates all the types of mutations"""
        return [getattr(MutationType, attr) for attr in filter(lambda a:
                a.startswith('MUT_'), dir(MutationType))]


class Mutations(object):
    """Keeps track of mutations for a given germline"""
    def __init__(self, germline, cdr3_nts):
        """Initializes the mutation statistics with a given germline & CDR3
        length."""
        self.region_stats = {}
        for region in ['all', 'CDR1', 'CDR2', 'CDR3', 'FR1', 'FR2', 'FR3']:
            self.region_stats[region] = self._create_count_record()
        self.pos_stats = {}
        self.cdr3_nts = cdr3_nts
        self.germline = germline[:309] + cdr3_nts + germline[309+len(cdr3_nts):]

    def _get_region(self, index):
        """Determines the gene region from an offset index"""
        if index <= 77:
            return 'FR1'
        elif index <= 113:
            return 'CDR1'
        elif index <= 164:
            return 'FR2'
        elif index <= 194:
            return 'CDR2'
        elif index <= 308:
            return 'FR3'
        elif index <= 308 + len(self.cdr3_nts):
            return 'CDR3'
        return 'FR4'

    def _create_count_record(self, int_count=False):
        """Creates a statistics record for a region or position.  If
        `int_count` is True, only the counts will be maintained, otherwise all
        values will be stored."""
        rec = {}
        for m in (MutationType.MUT_SYN, MutationType.MUT_CONS,
                  MutationType.MUT_UNCONS):
            if int_count:
                rec[MutationType.get_readable(m)] = 0
            else:
                rec[MutationType.get_readable(m)] = {}
        return rec

    def _add_region_stat(self, i, seq):
        """Adds mutations from `seq` at a given position to the region stats"""
        mtype = self._get_mut_type(seq, i)
        if mtype not in (MutationType.MUT_NONE, MutationType.MUT_UNK):
            mtype = MutationType.get_readable(mtype)
            region = self._get_region(i)
            if region not in self.region_stats:
                self.region_stats[region] = self._create_count_record()

            mutation = (i, self.germline[i], seq[i],
                        self._get_aa_at(self.germline, i),
                        self._get_aa_at(seq, i))
            if mutation not in self.region_stats[region][mtype]:
                self.region_stats[region][mtype][mutation] = 0
            self.region_stats[region][mtype][mutation] += 1

            if mutation not in self.region_stats['all'][mtype]:
                self.region_stats['all'][mtype][mutation] = 0
            self.region_stats['all'][mtype][mutation] += 1

    def _add_pos_stat(self, i, mtype, seq):
        """Adds mutations from `seq` at a given position to the position
        stats"""
        mtype = self._get_mut_type(seq, i)
        if mtype not in (MutationType.MUT_NONE, MutationType.MUT_UNK):
            if i not in self.pos_stats:
                self.pos_stats[i] = self._create_count_record(True)
            self.pos_stats[i][MutationType.get_readable(mtype)] += 1

    def _get_aa_at(self, seq, i):
        """Gets the amino acid that is partially encoded by position `i`"""
        aa_off = i - i % 3
        return lookups.aa_from_codon(seq[aa_off:aa_off + 3])

    def _get_mut_type(self, seq, i):
        """Determines the mutation type of a sequence at a position"""
        if self.germline[i] != seq[i]:
            grm_aa = self._get_aa_at(self.germline, i)
            seq_aa = self._get_aa_at(seq, i)

            if grm_aa is None or seq_aa is None:
                return MutationType.MUT_UNK
            elif grm_aa != seq_aa:
                if lookups.are_conserved_aas(grm_aa, seq_aa):
                    return MutationType.MUT_CONS
                return MutationType.MUT_UNCONS
            else:
                return MutationType.MUT_SYN
        else:
            return MutationType.MUT_NONE

    def _get_replacement(self, i, to):
        return self.germline[:i] + to + self.germline[i+1:]

    def add_sequence(self, seq):
        """Calculates all mutation information for a sequence"""
        mut_str = ''
        for i in range(0, len(seq)):
            mut = self._get_mut_type(seq, i)
            mut_str += MutationType.get_symbol(mut)
            self._add_region_stat(i, seq)
            self._add_pos_stat(i, mut, seq)
        return mut_str

    def get_aggregate(self):
        """Aggregates all mutation information from added sequences"""
        region_stats = {}
        for region, stats in self.region_stats.iteritems():
            region_stats[region] = {
                'counts': {
                    'unique': self._create_count_record(True),
                    'total': self._create_count_record(True)
                },
                'mutations': {}
            }

            for mut_type, mutations in stats.iteritems():
                region_stats[region]['counts']['total'][mut_type] = \
                    sum(mutations.values())
                region_stats[region]['counts']['unique'][mut_type] = \
                    len(mutations)

                region_stats[region]['mutations'][mut_type] = []
                for mutation, count in mutations.iteritems():
                    region_stats[region]['mutations'][mut_type].append({
                        'count': count,
                        'position': mutation[0],
                        'from': mutation[1],
                        'to': mutation[2],
                        'aa_from': mutation[3],
                        'aa_to': mutation[4],
                    })

        return region_stats, self.pos_stats
