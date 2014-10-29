from sldb.common.mutations import *

def _same_dict(d1, d2):
    return cmp(d1, d2) == 0

def test_germline_comparison():
    germline = '''
    CAGGTTCAGCTGGTGCAGTCTGGAGCT---GAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTCTCCTGCAAGGCT
    TCTGGTTACACCTTT------------ACCAGCTATGGTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTT
    GAGTGGATGGGATGGATCAGCGCTTAC------AATGGTAACACAAACTATGCACAGAAGCTCCAG---GGCAGA
    GTCACCATGACCACAGACACATCCACGAGCACAGCCTACATGGAGCTGAGGAGCCTGAGATCTGACGACACGGCC
    GTGTATTACNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNGGCCAGGGCACCCTGGTCACCGTC
    TCCTCAG
    '''.replace('\n', '').strip()

    seq = '''
    NNNNNNNNNNNNNNNNNNNNNNNNNNN---NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNCTGCAAGGCT
    TCTGGAGGCACCTTC------------AGCAGCTATGCTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTT
    GAGTGGATGGGATGGATCAGCGCTTAC------AATGGTAACACAAACTATGCACAGAAGCTCCAG---GGCAGA
    GTCACCATGACCACAGACACATCCACGAGCACAGCCTACATGGAGCTGAGGAGCCTGAGATCTGACGACACGGCC
    GTGTATTACTGTGCGTGAGCAGCAGCTGGTACGAATGGGTTCGACCCCTGGGGCCAGGGAACCCTGGTCACCGTC
    TCCTCAG
    '''.replace('\n', '').strip()

    cdr3_len = 42

    mut = Mutations(germline, cdr3_len)
    mut.add_sequence(seq)

    regions, poss = mut.get_aggregate()

    assert _same_dict(poss, 
        {  
            116: {  
                'synonymous':0,
                'nonconservative':0,
                'conservative':1
            },
            84: {  
                'synonymous':0,
                'nonconservative':1,
                'conservative':0
            },
            85: {  
                'synonymous':0,
                'nonconservative':1,
                'conservative':0
            },
            86: {  
                'synonymous':0,
                'nonconservative':1,
                'conservative':0
            },
            375: {  
                'synonymous':0,
                'nonconservative':1,
                'conservative':0
            }
        })
