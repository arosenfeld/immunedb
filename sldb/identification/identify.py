import argparse
from Bio import SeqIO

from sldb.identification.vdj_sequence import VDJSequence
from sldb.identification.v_ties import get_v_ties


def _run_sample(fh):
    lengths_sum = 0
    mutations_sum = 0
    vdjs = []

    no_result = 0
    for i, record in enumerate(SeqIO.parse(fh, 'fasta')):
        if i > 0 and i % 1000 == 0:
            print i
        vdj = VDJSequence(record.description, record.seq)
        if vdj.j_gene is not None and vdj.v_gene is not None:
            lengths_sum += vdj.v_anchor_pos
            mutations_sum += vdj.mutation_fraction
            vdjs.append(vdj)
        else:
            no_result += 1

    avg_len = lengths_sum / float(len(vdjs))
    avg_mutation_frac = mutations_sum / float(len(vdjs))
    print 'len={}, mut={}'.format(avg_len, avg_mutation_frac)
    v_ties = get_v_ties(avg_len, avg_mutation_frac)

    for vdj in vdjs:
        new_vs = set([])
        for v in vdj.v_gene:
            if v in v_ties:
                ties = v_ties[v].split('|')
                new_vs = new_vs.union(set(ties))
            else:
                new_vs.add(v)
        old_vs = vdj.v_gene[:]
        vdj.v_gene = list(new_vs)
        print vdj.id, vdj.v_gene, old_vs
        print vdj.germline
        print vdj.sequence

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta_file')
    args = parser.parse_args()

    with open(args.fasta_file) as fh:
        _run_sample(fh)
