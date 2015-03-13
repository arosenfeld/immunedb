import json

from sldb.common.models import Sequence
from sldb.identification.v_genes import VGene
import sldb.util.lookups as lookups


class ContextualMutations(object):
    def __init__(self):
        self._seen = {}
        self.region_muts = {}
        self.position_muts = {}

    def _get_region(self, index, cdr3_num_nts):
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
        elif index <= 308 + cdr3_num_nts:
            return 'CDR3'
        return 'FR4'

    def add_mutation(self, seq_replaced, cdr3_num_nts, mtype, mutation,
                     copy_number):
        pos, from_nt, to_nt, from_aa, to_aa = mutation
        uniq = '_'.join(map(str, (pos, from_nt, to_nt)))
        region = self._get_region(pos, cdr3_num_nts)
        if region not in self.region_muts:
            self.region_muts[region] = {
                'total': {},
                'unique': {}
            }
        # If it's a new mutation, setup the dictionaries
        if uniq not in self.region_muts[region]['total']:
            self.region_muts[region]['total'][uniq] = {}
            self.region_muts[region]['unique'][uniq] = {}

        # If the mutation type has not been seen for the mutation, setup the
        # dictionaries
        if mtype not in self.region_muts[region]['total'][uniq]:
            self.region_muts[region]['total'][uniq][mtype] = {}
            self.region_muts[region]['unique'][uniq][mtype] = {}

        # If the specific amino acid change has not been seen, setup the
        # counters
        aa_key = '_'.join((from_aa, to_aa))
        if aa_key not in self.region_muts[region]['total'][uniq][mtype]:
            self.region_muts[region]['total'][uniq][mtype][aa_key] = 0
            self.region_muts[region]['unique'][uniq][mtype][aa_key] = 0

        # Increment the total count for the mutation
        self.region_muts[region]['total'][uniq][mtype][aa_key] = copy_number

        # If this is the first time the sequence has been seen, add it to the
        # unique count.
        if uniq not in self._seen:
            self._seen[uniq] = set([])
        if seq_replaced not in self._seen[uniq]:
            self.region_muts[region]['unique'][uniq][mtype][aa_key] = 1
            self._seen[uniq].add(seq_replaced)
            if pos not in self.position_muts:
                self.position_muts[pos] = {}
            if mtype not in self.position_muts[pos]:
                self.position_muts[pos][mtype] = 0
            self.position_muts[pos][mtype] += 1


class CloneMutations(object):
    def __init__(self, session, clone):
        self._clone = clone
        self._session = session

        self._germline = self._clone.group.germline
        self._germline = ''.join([
            self._germline[:VGene.CDR3_OFFSET],
            clone.cdr3_nt,
            self._germline[VGene.CDR3_OFFSET + clone.cdr3_num_nts:]
        ])

    def _get_aa_at(self, seq, i):
        aa_off = i - i % 3
        return lookups.aa_from_codon(seq[aa_off:aa_off + 3])

    def _get_mutation(self, seq, i):
        if (self._germline[i] != seq[i]
                and self._germline[i] != 'N'
                and seq[i] != 'N'):
            grm_aa = self._get_aa_at(self._germline, i)
            seq_aa = self._get_aa_at(seq, i)

            if grm_aa is None or seq_aa is None:
                return 'unknown'
            elif grm_aa != seq_aa:
                if lookups.are_conserved_aas(grm_aa, seq_aa):
                    return 'conservative'
                return 'nonconservative'
            else:
                return 'synonymous'

        return None

    def calculate(self, commit_seqs=False, limit_samples=None):
        sample_mutations = {}
        clone_mutations = ContextualMutations()

        seqs = self._session.query(Sequence).filter(
            Sequence.clone == self._clone)
        if limit_samples is not None:
            seqs = seqs.filter(Sequence.sample_id.in_(limit_samples))

        for seq in seqs:
            seq_mutations = []
            for i in range(0, len(seq.sequence)):
                if seq.sample_id not in sample_mutations:
                    sample_mutations[seq.sample_id] = ContextualMutations()
                mtype = self._get_mutation(seq.sequence, i)
                if mtype is None:
                    continue
                mutation = (i, self._germline[i], seq.sequence[i],
                            self._get_aa_at(self._germline, i) or '?',
                            self._get_aa_at(seq.sequence, i) or '?')
                seq_mutations.append(mutation)

                sample_mutations[seq.sample_id].add_mutation(
                    seq.sequence_replaced, self._clone.cdr3_num_nts,
                    mtype, mutation, seq.copy_number)

                clone_mutations.add_mutation(
                        seq.sequence_replaced, self._clone.cdr3_num_nts, mtype,
                        mutation, seq.copy_number)
            if commit_seqs:
                seq.mutations_from_clone = json.dumps(seq_mutations)

        if commit_seqs:
            self._session.commit()
        return sample_mutations, clone_mutations
