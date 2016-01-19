import os

from Bio.Seq import Seq

from sldb.identification.vdj_sequence import VDJSequence
from sldb.identification.v_genes import VGermlines
from sldb.identification.j_genes import JGermlines
import sldb.util.lookups as lookups


def test_vdj():
    '''Tests VDJ identification on one sequence'''
    v_germs = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'imgt_human_v.fasta')
    j_germs = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'imgt_human_j.fasta')
    seq_id = 'M01651:98:000000000-A7TTB:1:1101:10239:17854'
    seq = ('TACGTACTGAGGAGACGGTGACCAGGGTTCCCTGGCCCCAGGGGTCGAACCAGTTGCCTC'
           'GGAGATTACTACTCCAAAAATCGTAATTCTCTCTCGCACAATAATAAATGGCCGTGTCCG'
           'CAGCGGTCGCCGAGCTCAGGTTCAGAAAGAACTGATTCTTGGACGTGTCCAGTGATATGG'
           'CGACTCGACCCTCGAGGGAGGGGTTGTAGTTGGTGCTCCCAGTGTAGTAGATGTGTCCAA'
           'TCCACTCCAGCCCCTTCCCTGGGGGCTGCCGGATCCAGTTCCAGTTGTAAGAACCGAAGG'
           'AGCCCCCAGAGACAGTGCAGGTGAGGGACAGGGTCTCCGAGGACTTCACCAGTCCTGGGC'
           'CCGACTCCTGCAGCTGCACCTGGGACAGGACCCCTGTGAACAGAGAGACCCACAGTGAGC'
           'CCTGGGATCAGAGGCACCTCCCATATCCCCATGTCTGGATCCCTGAGACACTCACATCTG'
           'GGAGCTGCCACCAGGAGAAGGAAGAACCACAGGTGTTTCAT')

    v_germlines = VGermlines(v_germs)
    j_germlines = JGermlines(j_germs, 31, 18, 12)
    vdj = VDJSequence(seq_id, seq, v_germlines, j_germlines, analyze=True)
    vdj.align_to_germline(vdj.v_length, vdj.mutation_fraction)
    assert set(vdj.v_gene) == set(['IGHV4-59*01', 'IGHV4-59*02',
                                   'IGHV4-59*08'])
    assert set(vdj.j_gene) == set([
        'IGHJ1*01', 'IGHJ4*01', 'IGHJ4*02', 'IGHJ4*03', 'IGHJ5*01', 'IGHJ5*02'
    ])
    assert vdj.cdr3 == ('TGTGCGAGAGAGAATTACGATTTTTGGAGTAGTAATCTCCGAGGCAACTGGT'
                        'TCGACCCCTGG')
    assert lookups.aas_from_nts(vdj.cdr3, '') == 'CARENYDFWSSNLRGNWFDPW'

    assert vdj.post_cdr3_length == 31
    assert vdj.v_match == 265
    assert vdj.j_match == 48
    assert vdj.post_cdr3_match == 31
    assert vdj.pre_cdr3_match == 281
    assert vdj.j_length == 48
    assert vdj.v_length == 293
    assert vdj.quality is None
    assert len(vdj.ids) == 1
    assert vdj.ids[0] == 'M01651:98:000000000-A7TTB:1:1101:10239:17854'
    assert vdj.removed_prefix_qual == ''
    assert vdj.germline == (
        'CAGGTGCAGCTGCAGGAGTCGGGCCCA---GGACTGGTGAAGCCTTCGGAGACCCTGTCCCTCACC'
        'TGCACTGTCTCTGGTGGCTCCNTC------------AGTAGTTACTACTGGAGCTGGATCCGGCAG'
        'CCCCCAGGGAAGGGACTGGAGTGGATTGGGTATATCTATTACAGT---------GGGAGCACCAACT'
        'ACAACCCCTCCCTCAAG---AGTCGAGTCACCATATCAGTAGACACGTCCAAGAACCAGTTCTCCCT'
        'GAAGCTGAGCTCTGTGACCGCNGCNGACACGGCCGTGTATTAC------------------------'
        '---------------------------------------GGCCANGGNACCCTGGTCACCGTCTCCT'
        'CAG'
    )
    assert vdj.pre_cdr3_length == 309
    assert vdj.sequence == (
        'CAGGTGCAGCTGCAGGAGTCGGGCCCA---GGACTGGTGAAGTCCTCGGAGACCCTGTCCCTCACCT'
        'GCACTGTCTCTGGGGGCTCCTTC------------GGTTCTTACAACTGGAACTGGATCCGGCAGCC'
        'CCCAGGGAAGGGGCTGGAGTGGATTGGACACATCTACTACACT---------GGGAGCACCAACTAC'
        'AACCCCTCCCTCGAG---GGTCGAGTCGCCATATCACTGGACACGTCCAAGAATCAGTTCTTTCTGA'
        'ACCTGAGCTCGGCGACCGCTGCGGACACGGCCATTTATTATTGTGCGAGAGAGAATTACGATTTTTG'
        'GAGTAGTAATCTCCGAGGCAACTGGTTCGACCCCTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCA'
        'G'
    )
    assert vdj.j_anchor_pos == 448
