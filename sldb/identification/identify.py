import argparse
from Bio import SeqIO

from sldb.identification.vdj_sequence import VDJSequence

def _run_sample(fh):
    mutation_fraction = []
    no_result = 0
    for i, record in enumerate(SeqIO.parse(fh, 'fasta')):
        if i > 0 and i % 1000 == 0:
            print i
            if i == 5000:
                break
        vdj = VDJSequence(record.seq)
        if vdj.j_gene is not None and vdj.v_gene is not None:
            mutation_fraction.append(vdj.mutation_fraction)
        else:
            no_result += 1
    print 'mutation', sum(mutation_fraction) / float(len(mutation_fraction))
    print 'results', len(mutation_fraction)
    print 'noresults', no_result

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta_file')
    args = parser.parse_args()

    with open(args.fasta_file) as fh:
        _run_sample(fh)
#'M01651:98:000000000-A7TTB:1:1101:23634:8217 1:N:0:27': # DC
#'M01651:98:000000000-A7TTB:1:1101:10142:17712 1:N:0:27': # YXC
#'M01651:98:000000000-A7TTB:1:1101:11042:6786 1:N:0:27': # Ties
