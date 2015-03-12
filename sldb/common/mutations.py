import json

from sldb.common.models import Sequence
from sldb.identification.v_genes import VGene
import sldb.util.lookups as lookups

class MutationCount(object):
    def __init__(self, session, clone, seqs=None):
        self._clone = clone
        self._session = session

        if seqs is None:
            self._seqs = session.query(Sequence).filter(
                Sequence.clone == clone).all()
        else:
            self._seqs = seqs

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

        for seq in self._seqs:
            seq_mutations = {}
            for i in range(0, len(seq.sequence)):
                mtype = self._get_mutation(seq.sequence, i)
                if mtype is not None:
                    # String for jsonifying
                    key = '_'.join(map(str, (i, self._germline[i],
                                             seq.sequence[i])))
                    if key not in seq_mutations:
                        seq_mutations[key] = {
                            'type': mtype,
                        }

                    if seq.sample_id not in sample_mutations:
                        sample_mutations[seq.sample_id]['mutations'] = {
                            'unique': {},
                            'total': {},
                            'mutations': {}
                        }
                    if key not in sample_mutations[seq.sample_id]:
                        sample_mutations[seq.sample_id]['mutations'][key] = {
                            'type': mtype,
                            'unique': 0,
                            'total': 0,
                        }
                if commit_seqs:
                    seq.mutations_from_clone = json.dumps(seq_mutations)

            # Update the sample mutations
            sample_mutations[seq.sample_id]['mutations'][key]['unique'] += 1
            sample_mutations[seq.sample_id]['mutations'][key]['total'] += (
                    seq.copy_number)
            # Update the sample counts
            sample_mutations[seq.sample_id]['unique'][mtype] += 1
            sample_mutations[seq.sample_id]['total'][mtype] += seq.copy_number

        if commit_seqs:
            self._session.commit()

        return sample_mutations
