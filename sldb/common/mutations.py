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

    def add_mutation(self, seq_replaced, cdr3_num_nts, mutation, from_aa,
                     intermediate_seq_aa, final_seq_aa, copy_number):
        pos, from_nt, to_nt, mtype = mutation
        region = self._get_region(pos, cdr3_num_nts)
        if region not in self.region_muts:
            self.region_muts[region] = {
                'counts': {
                    'total': {},
                    'unique': {}
                },
                'mutations': {}
            }

        # If it's a new mutation, setup the dictionaries
        if mtype not in self.region_muts[region]['counts']['total']:
            self.region_muts[region]['counts']['total'][mtype] = 0
            self.region_muts[region]['counts']['unique'][mtype] = 0
            self.region_muts[region]['mutations'][mtype] = {}

        if mutation not in self.region_muts[region]['mutations'][mtype]:
            self.region_muts[region]['mutations'][mtype][mutation] = {
                'pos': pos,
                'from_nt': from_nt,
                'from_aa': from_aa,
                'to_nt': to_nt,
                'to_aas': [],

                'unique': 0,
                'total': 0,
                'intermediate_aa': intermediate_seq_aa,
            }

        mut_dict = self.region_muts[region]['mutations'][mtype][mutation]
        mut_dict['total'] += copy_number
        if final_seq_aa not in mut_dict['to_aas']:
            mut_dict['to_aas'].append(final_seq_aa)

        self.region_muts[region]['counts']['total'][mtype] += copy_number

        if mutation not in self._seen:
            self._seen[mutation] = set([])
        if seq_replaced not in self._seen[mutation]:
            self.region_muts[region]['counts']['unique'][mtype] += 1
            mut_dict['unique'] += 1
            self._seen[mutation].add(seq_replaced)

    def get_all(self):
        final_regions = {}
        for region, stats in self.region_muts.iteritems():
            final_regions[region] = {
                'counts': stats['counts'],
                'mutations': {}
            }

            for mtype, muts in stats['mutations'].iteritems():
                final_regions[region]['mutations'][mtype] = muts.values()

        return {
            'regions': final_regions,
            'positions': self.position_muts
        }


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

    def _get_codon_at(self, seq, i):
        aa_off = i - i % 3
        return seq[aa_off:aa_off + 3]

    def _get_aa_at(self, seq, i):
        return lookups.aa_from_codon(self._get_codon_at(seq, i))

    def _get_mutation(self, seq, i):
        if (self._germline[i] != seq[i]
                and self._germline[i] != 'N'
                and seq[i] != 'N'):
            grm_aa = self._get_aa_at(self._germline, i)
            # Simulate this mutation alone
            off = i % 3
            grm_codon = self._get_codon_at(self._germline, i)
            seq_aa = lookups.aa_from_codon(
                grm_codon[:off] + seq[i] + grm_codon[off+1:])

            if grm_aa is None or seq_aa is None:
                return 'unknown', seq_aa
            elif grm_aa != seq_aa:
                if lookups.are_conserved_aas(grm_aa, seq_aa):
                    return 'conservative', seq_aa
                return 'nonconservative', seq_aa
            else:
                return 'synonymous', seq_aa

        return None, None

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
                mtype, intermediate_seq_aa = self._get_mutation(
                    seq.sequence, i)
                if mtype is None:
                    continue

                mutation = (i, self._germline[i], seq.sequence[i], mtype)
                from_aa = self._get_aa_at(self._germline, i)
                seq_mutations.append(mutation)

                sample_mutations[seq.sample_id].add_mutation(
                    seq.sequence_replaced, self._clone.cdr3_num_nts, mutation,
                    from_aa, intermediate_seq_aa,
                    self._get_aa_at(seq.sequence, i), seq.copy_number)

                clone_mutations.add_mutation(
                    seq.sequence_replaced, self._clone.cdr3_num_nts, mutation,
                    from_aa, intermediate_seq_aa,
                    self._get_aa_at(seq.sequence, i), seq.copy_number)
            if commit_seqs:
                seq.mutations_from_clone = json.dumps(seq_mutations)

        if commit_seqs:
            self._session.commit()
        return sample_mutations, clone_mutations

def create_threshold_mutations(all_muts, total_seqs, threshold_type,
                               threshold_value):
    pas
