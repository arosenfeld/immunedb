import argparse
import pickle

from Bio.Seq import Seq

from sldb.identification.vdj_sequence import VDJSequence
from sldb.identification.v_genes import VGermlines
from sldb.identification.j_genes import JGermlines

def compare_dicts(regression, new):
    assert set(regression.keys()) == set(new.keys())
    for k in regression.keys():
        assert regression[k] == new[k], (
            'Mismatch in sequence {} for key {}, reg={}, new={}'.format(
                new['id'], k, regression[k], new[k]))

parser = argparse.ArgumentParser()
parser.add_argument('regression_file')
parser.add_argument('v_germline_file')
parser.add_argument('j_germline_file')
parser.add_argument('--partial-reads', default=False, action='store_true')

args = parser.parse_args()

fields = [
    'id',
    'j_gene',
    'v_gene',
    'j_anchor_pos',
    'v_anchor_pos',
    'cdr3',
    'sequence',
    'sequence_filled',
    'germline',
    'mutation_fraction',
    'num_gaps',
    'pad_length',
    'in_frame',
    'stop',
    'functional',
    'j_length',
    'j_match',
    'v_length',
    'v_match',
    'pre_cdr3_length',
    'pre_cdr3_match',
    'post_cdr3_length',
    'post_cdr3_match',
]

v_germlines = VGermlines(args.v_germline_file)
j_germlines = JGermlines(args.j_germline_file, upstream_of_cdr3=31,
                         anchor_len=18, min_anchor_len=12)
with open(args.regression_file) as reg_file:
    regression = pickle.load(reg_file)

for i, record in enumerate(regression):
    vdj = VDJSequence(record['id'], record['original_seq'], not
                      args.partial_reads, v_germlines, j_germlines)
    new = {'id': vdj.id, 'v_gene': None, 'j_gene': None}
    if vdj.v_gene is not None and vdj.j_gene is not None:
        vdj.align_to_germline()
        if vdj.v_gene is not None and vdj.j_gene is not None:
            new = {field: getattr(vdj, field) for field in fields}
    new['original_seq'] = record['original_seq']

    compare_dicts(record, new)
    print '{} / {}'.format(i + 1, len(regression))
