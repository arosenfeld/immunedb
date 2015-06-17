import argparse
import pickle

from Bio.Seq import Seq

from sldb.identification.vdj_sequence import AlignmentException, VDJSequence
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
    new = {'id': record['id'], 'v_gene': None, 'j_gene': None}
    try:
        vdj = VDJSequence(record['id'], record['original_seq'], v_germlines,
                j_germlines)
        vdj.align_to_germline()
        new = {field: getattr(vdj, field) for field in fields}
    except AlignmentException:
        pass
    new['original_seq'] = record['original_seq']

    compare_dicts(record, new)
    if i % 1000 == 0:
        print '{} / {}'.format(i + 1, len(regression))
