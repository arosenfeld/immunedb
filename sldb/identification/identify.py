import json
import os
import re
from Bio import SeqIO

from sldb.identification.vdj_sequence import VDJSequence
from sldb.identification.v_ties import get_v_ties

def _identify_reads(fn):
    lengths_sum = 0
    mutations_sum = 0
    vdjs = []

    no_result = 0
    print 'Starting {}'.format(fn)
    for i, record in enumerate(SeqIO.parse(fn, 'fasta')):
        if i > 0 and i % 1000 == 0:
            print i
        vdj = VDJSequence(record.description, record.seq)
        if vdj.j_gene is not None and vdj.v_gene is not None:
            lengths_sum += vdj.v_anchor_pos
            mutations_sum += vdj.mutation_fraction
            vdjs.append(vdj)
        else:
            no_result += 1
    cnt = i

    avg_len = lengths_sum / float(len(vdjs))
    avg_mutation_frac = mutations_sum / float(len(vdjs))
    v_ties = get_v_ties(avg_len, avg_mutation_frac)

    v_tie_cnt = 0
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
        if len(vdj.v_gene) > 1:
            v_tie_cnt += 1
    print 'total_seqs={}'.format(cnt)
    print 'vties={}'.format(v_tie_cnt)
    print 'len={}'.format(avg_len)
    print 'mut={}'.format(avg_mutation_frac)
    print 'noresults={}'.format(no_result)


def run_identify(args):
    with open('{}/metadata.json'.format(args.base_dir)) as fh:
        metadata = json.load(fh)

        names = set([])
        for fn in os.listdir('{}/processed'.format(args.base_dir)):
            name = fn.split('.')[0].rsplit('_', 2)[0]
            names.add(name)

        for name in names:
            base = '{}/presto/{}'.format(args.base_dir, name)
            r1 = '{}_R1_001.sync_assemble-fail.fasta'.format(base)
            r2 = '{}_R2_001.sync_assemble-fail.fasta'.format(base)
            join = '{}_R1_001.sync_assemble-pass.fasta'.format(base)
            if not os.path.isfile('{}.log'.format(base)):
                print 'Skipping {} since no presto log exists.'.format(name)
                continue
            _identify_reads(join)
