import sys
from Bio import SeqIO
from Bio.Seq import Seq

from sldb.identification.vdj_sequence import VDJSequence
from sldb.identification.v_genes import VGene, VGermlines
import sldb.identification.v_ties as v_ties
import sldb.util.lookups as lookups


if __name__ == '__main__':
    v_germlines = VGermlines(sys.argv[1])
    for record in SeqIO.parse(sys.argv[2], 'fasta'):
        print record.description
        vdj = VDJSequence(record.description, record.seq, 'presto' in
                          record.description, v_germlines)
        vdj.align_to_germline(vdj.v_length, vdj.mutation_fraction)
        if vdj.j_gene is not None and vdj.v_gene is not None:
            print 'v_gene     :', vdj.v_gene
            print 'j_gene     :', vdj.j_gene
            print 'cdr3 len   :', len(vdj.cdr3)
            print 'cdr3       :', vdj.cdr3
            print 'cdr3_aa    :', lookups.aas_from_nts(vdj.cdr3, '')
            print 'germ       :', vdj.germline
            print 'seq        :', vdj.sequence
            print 'v_length   :', vdj.v_length
            print 'v_match    :', vdj.v_match
            print 'pre_length :', vdj.pre_cdr3_length
            print 'pre_match  :', vdj.pre_cdr3_match
            print ''
        else:
            print 'No result'
