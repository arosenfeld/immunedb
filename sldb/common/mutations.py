import json

from sldb.common.models import Sequence
from sldb.identification.v_genes import VGene
import sldb.util.lookups as lookups

class ContextualMutations(object):
    def __init__(self):
        self._seen = set([])
        self._mutations = {
            'total': {},
            'unique': {},
        }

    def add_mutation(self, seq_replaced, mtype, mutation, copy_number):
        pos, from_nt, to_nt, from_aa, to_aa = mutation
        uniq = (pos, from_nt, to_nt)
        # If it's a new mutation, setup the dictionaries
        if uniq not in self['total']:
            self._mutations['total'][uniq] = {}
            self._mutations['unique'][uniq] = {}

        # If the mutation type has not been seen for the mutation, setup the
        # dictionaries
        if mtype not in record['total'][uniq]:
            record['total'][uniq][mtype] = {}
            record['unique'][uniq][mtype] = {}

        # If the specific amino acid change has not been seen, setup the
        # counters
        if (from_aa, to_aa) not in record['total'][uniq][mtype]:
            record['total'][uniq][mtype][(from_aa, to_aa)] = 0
            record['unique'][uniq][mtype][(from_aa, to_aa)] = 0

        # Increment the total count for the mutation
        record['total'][uniq][mtype][(from_aa, to_aa)] = copy_number

        # If this is the first time the sequence has been seen, add it to the
        # unique count.
        if seq_replaced not in self._seen:
            record['unique'][uniq][mtype][(from_aa, to_aa)] = 1
            self._seen.add(seq_replaced)

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
                    return 'conserved'
                return 'unconserved'
            else:
                return 'synonymous'

        return None

    def calculate(self, commit_seqs=False):
        sample_mutations = {}
        total_mutations = ContextualMutations()

        for seq in session.query(Sequence).filter(
                Sequence.clone == clone):
            seq_mutations = []
            for i in range(0, len(seq.sequence)):
                mtype = self._get_mutation(seq.sequence, i)
                if mtype is None:
                    continue
                seq_mutations.add((i, self._germline[i], seq.sequence[i],
                                   from_aa, to_aa))

                mutation = (i, self._germline[i], seq.sequence[i],
                            self._get_aa_at(self._germline, i) or '?',
                            self._get_aa_at(seq.sequence, i) or '?')
                if seq.sample_id not in sample_mutations:
                    sample_mutations[seq.sample_id] = ContextualMutations()
                sample_mutations[seq.sample_id].add_mutation(
                        seq.sequence_replaced, mtype, mutation, seq.copy_number)

                total_mutations.add_mutation(seq.sequence_replaced, mtype,
                                             mutation, seq.copy_number)
            if commit_seqs:
                seq.mutations_from_clone = json.dumps(seq_mutations)

        if commit_seqs:
            self._session.commit()

        return sample_mutations, all_mutations
