from Bio.Seq import Seq

j_anchors = {
    'IGHJ1|4|5': Seq('TGGTCACCGTCTCCTCAG'),
    'IGHJ2':     Seq('TGGTCACTGTCTCCTCAG'),
    'IGHJ3':     Seq('TGGTCACCGTCTCTTCAG'),
    'IGHJ6':     Seq('CGGTCACCGTCTCCTCAG'),
}

j_anchors = {
    'IGHJ1':
    Seq('-----------GCTGAATACTTCCAGCACTGGGGCCAGGGCACCCTGGTCACCGTCTCCTCAG'),
    'IGHJ4':
    Seq('---------------ACTACTTTGACTACTGGGGCCAAGGAACCCTGGTCACCGTCTCCTCAG'),
    'IGHJ5':
    Seq('------------ACAACTGGTTCGACTCCTGGGGCCAAGGAACCCTGGTCACCGTCTCCTCAG'),
    'IGHJ2':
    Seq('----------CTACTGGTACTTCGATCTCTGGGGCCGTGGCACCCTGGTCACTGTCTCCTCAG'),
    'IGHJ3':
    Seq('-------------TGATGCTTTTGATGTCTGGGGCCAAGGGACAATGGTCACCGTCTCTTCAG'),
    'IGHJ6':
    Seq('ATTACTACTACTACTACGGTATGGACGTCTGGGGGCAAGGGACCACGGTCACCGTCTCCTCAG'),
}


v_regex = [
    'GAC.{15}TGT',
    'GAC.{15}TG',
]


dc_final_aas = ['YY', 'YC', 'YH']


germline_anchor = 'TGT'


def all_j_anchors():
    max_size = max(map(len, j_anchors.values()))
    for trim in range(0, max_size, 3):
        for j, f in j_anchors.iteritems():
            new_f = f[0:len(f) - trim]
            if len(new_f) >= 12:
                yield new_f, j
