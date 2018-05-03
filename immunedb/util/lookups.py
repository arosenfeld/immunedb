import itertools


def aas_from_nts(nts, replace_unknowns='X'):
    """Returns the amino acids for a given nuceotide string"""
    aas = ''
    for i in range(0, len(nts), 3):
        codon = str(nts[i:i+3])
        if len(codon) < 3:
            break
        aas += aa_from_codon(codon, replace_unknowns)
    return aas


def aa_from_codon(codon, ret_on_unknown=None):
    """Looks up an amino acid from a codon, or returns None"""
    if codon.upper() in _aa_lookup:
        return _aa_lookup[codon.upper()]
    return ret_on_unknown


def are_conserved_aas(aa1, aa2):
    """Determines if two amino acids are conserved based on the paper
    'Structural Determinants in the Sequences of Immunoglobulin Variable
    Domain'"""
    groups = [
        'WFILMVC',
        'PYTHSAG',
        'KEQNRD'
    ]
    for g in groups:
        if aa1.upper() in g and aa2.upper() in g:
            return True
    return False


def aa_to_all_nts(aas):
    nts = []
    for aa in aas:
        nts.append(_nt_lookup[aa])
    return [''.join(arr) for arr in itertools.product(*nts)]


def has_stop(seq):
    for i in range(0, len(seq), 3):
        if seq[i:i+3] in ['TAG', 'TAA', 'TGA']:
            return True
    return False


_aa_lookup = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
              'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y',
              'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*',
              'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
              'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H',
              'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R',
              'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
              'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
              'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S',
              'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V',
              'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A',
              'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
              'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}


_nt_lookup = {
    'A': ['GCA', 'GCC', 'GCG', 'GCT'],
    'C': ['TGT', 'TGC'],
    'E': ['GAG', 'GAA'],
    'D': ['GAT', 'GAC'],
    'G': ['GGT', 'GGG', 'GGA', 'GGC'],
    'F': ['TTT', 'TTC'],
    'I': ['ATC', 'ATA', 'ATT'],
    'H': ['CAT', 'CAC'],
    'K': ['AAG', 'AAA'],
    '*': ['TAG', 'TAA', 'TGA'],
    'M': ['ATG'],
    'L': ['CTT', 'CTG', 'CTA', 'CTC', 'TTA', 'TTG'],
    'N': ['AAC', 'AAT'],
    'Q': ['CAA', 'CAG'],
    'P': ['CCT', 'CCG', 'CCA', 'CCC'],
    'S': ['AGC', 'AGT', 'TCT', 'TCG', 'TCC', 'TCA'],
    'R': ['AGG', 'AGA', 'CGA', 'CGG', 'CGT', 'CGC'],
    'T': ['ACA', 'ACG', 'ACT', 'ACC'],
    'W': ['TGG'],
    'V': ['GTA', 'GTC', 'GTG', 'GTT'],
    'Y': ['TAT', 'TAC']
}
