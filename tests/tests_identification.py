import os

from Bio.Seq import Seq

from sldb.identification.vdj_sequence import VDJSequence
from sldb.identification.v_genes import VGermlines
import sldb.util.lookups as lookups

def test_full_seq():
    germs = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         'imgt_human_alleles.fasta')
    seq_id = 'M01651:98:000000000-A7TTB:1:1101:10239:17854'
    seq = Seq(
        'TACGTACTGAGGAGACGGTGACCAGGGTTCCCTGGCCCCAGGGGTCGAACCAGTTGCCTC'
        'GGAGATTACTACTCCAAAAATCGTAATTCTCTCTCGCACAATAATAAATGGCCGTGTCCG'
        'CAGCGGTCGCCGAGCTCAGGTTCAGAAAGAACTGATTCTTGGACGTGTCCAGTGATATGG'
        'CGACTCGACCCTCGAGGGAGGGGTTGTAGTTGGTGCTCCCAGTGTAGTAGATGTGTCCAA'
        'TCCACTCCAGCCCCTTCCCTGGGGGCTGCCGGATCCAGTTCCAGTTGTAAGAACCGAAGG'
        'AGCCCCCAGAGACAGTGCAGGTGAGGGACAGGGTCTCCGAGGACTTCACCAGTCCTGGGC'
        'CCGACTCCTGCAGCTGCACCTGGGACAGGACCCCTGTGAACAGAGAGACCCACAGTGAGC'
        'CCTGGGATCAGAGGCACCTCCCATATCCCCATGTCTGGATCCCTGAGACACTCACATCTG'
        'GGAGCTGCCACCAGGAGAAGGAAGAACCACAGGTGTTTCAT')

    v_germlines = VGermlines(germs)
    vdj = VDJSequence(seq_id, seq, True, v_germlines)
    vdj.align_to_germline(vdj.v_length, vdj.mutation_fraction)
    assert set(vdj.v_gene) == set(['IGHV4-59*01', 'IGHV4-59*02',
                                   'IGHV4-59*08'])
    assert set(vdj.j_gene) == set(['IGHJ1', 'IGHJ4', 'IGHJ5'])
    assert vdj.cdr3 == ('TGTGCGAGAGAGAATTACGATTTTTGGAGTAGTAATCTCCGAGGCAACTGGT'
                        'TCGACCCCTGG')
    assert lookups.aas_from_nts(vdj.cdr3, '') == 'CARENYDFWSSNLRGNWFDPW'

    assert vdj.v_length == 293
    assert vdj.v_match == 264
    assert vdj.pre_cdr3_length == 282
    assert vdj.pre_cdr3_match == 251
'''
v_length   : 293
v_match    : 264
pre_length : 282
pre_match  : 251
'''
