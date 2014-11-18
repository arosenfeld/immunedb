import itertools


def aas_from_nts(nts, ret_on_wrong=None):
    """Returns the amino acids for a given nuceotide string"""
    aas = ''
    for i in range(0, len(nts), 3):
        aas += aa_from_codon(str(nts[i:i+3]), ret_on_wrong)
    return aas


def aa_from_codon(codon, ret_on_wrong=None):
    """Looks up an amino acid from a codon, or returns None"""
    if codon.upper() in _aa_lookup:
        return _aa_lookup[codon.upper()]
    return ret_on_wrong


def aa_to_all_nts(aas):
    nts = []
    for aa in aas:
        nts.append(_nt_lookup[aa])
    return map(lambda arr: ''.join(arr), itertools.product(*nts))


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


_nt_lookup = {'A': ['GCA', 'GCC', 'GCG', 'GCT'], 'C': ['TGT', 'TGC'], 'E':
['GAG', 'GAA'], 'D': ['GAT', 'GAC'], 'G': ['GGT', 'GGG', 'GGA', 'GGC'], 'F':
['TTT', 'TTC'], 'I': ['ATC', 'ATA', 'ATT'], 'H': ['CAT', 'CAC'], 'K': ['AAG',
'AAA'], '*': ['TAG', 'TAA', 'TGA'], 'M': ['ATG'], 'L': ['CTT', 'CTG', 'CTA',
'CTC', 'TTA', 'TTG'], 'N': ['AAC', 'AAT'], 'Q': ['CAA', 'CAG'], 'P': ['CCT',
'CCG', 'CCA', 'CCC'], 'S': ['AGC', 'AGT', 'TCT', 'TCG', 'TCC', 'TCA'], 'R':
['AGG', 'AGA', 'CGA', 'CGG', 'CGT', 'CGC'], 'T': ['ACA', 'ACG', 'ACT', 'ACC'],
'W': ['TGG'], 'V': ['GTA', 'GTC', 'GTG', 'GTT'], 'Y': ['TAT', 'TAC']}


"""The ties between V's"""
v_gene_ties = {
    'IGHV1-18': 'IGHV1-18',
    'IGHV1-2': 'IGHV1-2',
    'IGHV1-24': 'IGHV1-24',
    'IGHV1-3': 'IGHV1-3',
    'IGHV1-45': 'IGHV1-45',
    'IGHV1-46': 'IGHV1-46',
    'IGHV1-58': 'IGHV1-58',
    'IGHV1-69': 'IGHV1-69',
    'IGHV1-8': 'IGHV1-8',
    'IGHV1-f': 'IGHV1-f',
    'IGHV2-26': 'IGHV2-26',
    'IGHV2-5': 'IGHV2-5',
    'IGHV2-70': 'IGHV2-70',
    'IGHV3-11': 'IGHV3-11|IGHV3-48',
    'IGHV3-13': 'IGHV3-13',
    'IGHV3-15': 'IGHV3-15',
    'IGHV3-20': 'IGHV3-20',
    'IGHV3-21': 'IGHV3-21',
    'IGHV3-23': 'IGHV3-23',
    'IGHV3-30': 'IGHV3-30|IGHV3-30-3|IGHV3-33',
    'IGHV3-30-3': 'IGHV3-30-3|IGHV3-30',
    'IGHV3-33': 'IGHV3-33|IGHV3-30',
    'IGHV3-43': 'IGHV3-43',
    'IGHV3-48': 'IGHV3-48|IGHV3-11',
    'IGHV3-49': 'IGHV3-49',
    'IGHV3-53': 'IGHV3-53|IGHV3-66',
    'IGHV3-64': 'IGHV3-64',
    'IGHV3-66': 'IGHV3-66|IGHV3-53',
    'IGHV3-7': 'IGHV3-7',
    'IGHV3-72': 'IGHV3-72',
    'IGHV3-73': 'IGHV3-73',
    'IGHV3-74': 'IGHV3-74',
    'IGHV3-9': 'IGHV3-9',
    'IGHV3-NL1': 'IGHV3-NL1',
    'IGHV3-d': 'IGHV3-d',
    'IGHV4-28': 'IGHV4-28',
    'IGHV4-30-2': 'IGHV4-30-2|IGHV4-30-4|IGHV4-39',
    'IGHV4-30-4': 'IGHV4-30-4|IGHV4-30-2',
    'IGHV4-31': 'IGHV4-31|IGHV4-59',
    'IGHV4-34': 'IGHV4-34',
    'IGHV4-39': 'IGHV4-39|IGHV4-59|IGHV4-b|IGHV4-30-2',
    'IGHV4-4': 'IGHV4-4|IGHV4-59|IGHV4-61',
    'IGHV4-59': 'IGHV4-59|IGHV4-61|IGHV4-31|IGHV4-39|IGHV4-4',
    'IGHV4-61': 'IGHV4-61|IGHV4-4|IGHV4-59',
    'IGHV4-b': 'IGHV4-b|IGHV4-39',
    'IGHV5-51': 'IGHV5-51',
    'IGHV5-a': 'IGHV5-a',
    'IGHV6-1': 'IGHV6-1',
    'IGHV7-4-1': 'IGHV7-4-1'
}
