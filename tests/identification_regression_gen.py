import argparse
import pickle

from Bio import SeqIO

from sldb.identification.vdj_sequence import AlignmentException, VDJSequence
from sldb.identification.v_genes import VGermlines
from sldb.identification.j_genes import JGermlines

parser = argparse.ArgumentParser()
parser.add_argument('fasta_file')
parser.add_argument('v_germline_file')
parser.add_argument('j_germline_file')
parser.add_argument('out_file')
parser.add_argument('--format', default='fastq')

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

regression = []
v_germlines = VGermlines(args.v_germline_file)
j_germlines = JGermlines(args.j_germline_file, upstream_of_cdr3=31,
                         anchor_len=18, min_anchor_len=12)

total = len(list(SeqIO.parse(args.fasta_file, args.format)))
for i, record in enumerate(SeqIO.parse(args.fasta_file, args.format)):
    vals = {'id': record.description, 'v_gene': None, 'j_gene': None}
    try:
        vdj = VDJSequence(record.description, str(record.seq), v_germlines,
                          j_germlines)
        vdj.align_to_germline()
        vals = {field: getattr(vdj, field) for field in fields}
    except AlignmentException:
        pass
    vals['original_seq'] = record.seq
    regression.append(vals)
    print '{} / {}'.format(i + 1, total)

with open(args.out_file, 'w+') as fh:
    pickle.dump(regression, fh)
