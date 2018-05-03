from collections import OrderedDict
import csv
import re

from immunedb.common.models import Sample, Sequence
from immunedb.util.log import logger


def _get_entries(seq, inference, gene_db):
    entry = OrderedDict()
    trim_start = re.search('[N-]*', seq.sequence).end()
    trim_end = re.search('[N-]*', seq.sequence[::-1]).end()
    trimmed_seq = seq.sequence[trim_start:len(seq.sequence) - trim_end]
    # Offsets of regions/segments
    v_gaps = trimmed_seq.count('-')
    trimmed_seq = trimmed_seq.replace('-', '')
    v_region = (1, len(trimmed_seq), 'V_region')
    v_segment = (1, 309 - v_gaps - trim_start, 'V_segment')
    cds = (
        v_segment[1] + 1,
        v_segment[1] + seq.cdr3_num_nts,
        'CDS'
    )
    j_segment = (cds[1] + 1, v_region[1], 'J_segment')

    # Key
    entry[('>Feature', seq.seq_id.split(' ')[0])] = []
    # V region
    entry[v_region] = []
    # V segment
    first_v = seq.v_gene.split('|')[0]
    entry[v_segment] = (
        ('gene', first_v),
        ('db_xref', '{}:{}'.format(gene_db, first_v)),
        ('inference', inference),
    )

    # J segment
    first_j = seq.j_gene.split('|')[0]
    entry[j_segment] = (
        ('gene', first_j),
        ('db_xref', '{}:{}'.format(gene_db, first_j)),
        ('inference', inference),
    )

    # CDS
    cds = ('<{}'.format(cds[0]), '>{}'.format(cds[1]), 'CDS')
    entry[cds] = (
        ('function', 'junction'),
        ('codon_start', 1),
        ('product', 'immunoglobulin heavy chain variable domain'),
        ('inference', inference)
    )

    return entry, trimmed_seq


def _write_sample(session, sample_id, gene_db, inference, header):
    for sample in session.query(Sample):
        seqs = session.query(
            Sequence.seq_id,
            Sequence.v_gene,
            Sequence.j_gene,
            Sequence.cdr3_num_nts,
            Sequence.sequence,
            Sequence.copy_number
        ).filter(
            Sequence.sample_id == sample_id
        ).filter(Sequence.stop == 0)
        with open('{}.tbl'.format(sample.name), 'w+') as gb_fh:
            writer = csv.writer(gb_fh, delimiter='\t')
            with open('{}.fsa'.format(sample.name), 'w+') as fasta_fh:
                for seq in seqs:
                    gb_entry, fasta_seq = _get_entries(
                        seq, inference, gene_db)
                    for entry, indented in gb_entry.items():
                        writer.writerow(entry)
                        if indented:
                            for indent in indented:
                                writer.writerow(('', '', '') + indent)
                    seq_header = header + ' [note=AIRR_READ_COUNT:{}]'.format(
                        seq.copy_number)
                    fasta_fh.write('>{}\n{}\n'.format(
                        seq.seq_id.split(' ')[0] + ' {}'.format(seq_header),
                        fasta_seq))


def write_genbank(session, args):
    args.inference = 'alignment:' + args.inference

    header = (
        '[organism={}] '
        '[moltype={}] '
        '[keywords=AIRR]').format(
        args.species, args.mol_type)

    for sample in session.query(Sample):
        logger.info('Exporting sample {}'.format(sample.name))
        _write_sample(session, sample.id, args.gene_db, args.inference,
                      header)
