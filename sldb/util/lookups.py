def aas_from_nts(nts):
    aas = ''
    for i in range(0, len(nts), 3):
        aas += aa_from_codon(nts[i:i+3])
    return aas

def aa_from_codon(codon):
    """Looks up an amino acid from a codon, or returns None"""
    aa_lookup = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S',
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

    if codon.upper() in aa_lookup:
        return aa_lookup[codon.upper()]
    return None


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
