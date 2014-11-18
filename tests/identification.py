import sys
from Bio import SeqIO
from sldb.identification.vdj_sequence import VDJSequence
import sldb.util.lookups as lookups


if __name__ == '__main__':
    for record in SeqIO.parse(sys.argv[1], 'fasta'):
        vdj = VDJSequence(record.description, record.seq)
        if vdj.j_gene is not None and vdj.v_gene is not None:
            if 'N' in vdj.sequence:
                print 'germ   :', vdj.germline
                print 'seq    :', vdj.sequence
                print 'seq_rep:', vdj.sequence_filled
                print 'cdr3   :', str(vdj.cdr3)
                print 'cdr3_aa:', lookups.aas_from_nts(vdj.cdr3, '')
                print 'v_gene :', vdj.v_gene
                print 'j_gene :', vdj.j_gene
                print ''
