import csv
import re

from immunedb.common.models import NoResult, Sample, Study, Subject
from immunedb.identification import (add_as_noresult, add_uniques,
                                     AlignmentException)
from immunedb.identification.vdj_sequence import VDJSequence
from immunedb.identification.j_genes import JGermlines
from immunedb.identification.v_genes import VGermlines
import immunedb.util.funcs as funcs
from immunedb.util.log import logger


DEFAULT_MAPPINGS = {
    'full_sequence': 'V-D-J-REGION',
    'seq_id': 'Sequence ID',
    'v_gene': 'V-GENE and allele',
    'j_gene': 'J-GENE and allele',
    'copy_number': None
}


def _collapse_seqs(session, sample, reader, columns):
    seqs = {}
    for record in reader:
        if record[columns.full_sequence] is None:
            NoResult(seq_id=record[columns.seq_id], sample_id=sample.id)
            continue

        record[columns.full_sequence] = record[columns.full_sequence].replace(
            '.', '').upper()
        if record[columns.full_sequence] not in seqs:
            seqs[record[columns.full_sequence]] = {
                'record': record,
                'seq_ids': []
            }
        seqs[record[columns.full_sequence]]['seq_ids'].append(
            record[columns.seq_id]
        )
        if columns.copy_number in record:
            for dup in record[columns.copy_number]:
                seqs[record[columns.full_sequence]]['seq_ids'].append(
                    '{}_DUP_{}'.format(columns.seq_id, dup)
                )

    return seqs.values()


def read_file(session, handle, sample, v_germlines, j_germlines, columns,
              remaps):
    seqs = _collapse_seqs(session, sample, csv.DictReader(handle,
                          delimiter='\t'), columns)

    aligned_seqs = {}
    missed = 0
    total = 0
    for total, seq in enumerate(seqs):
        if total > 0 and total % 1000 == 0:
            logger.info('Finished {}'.format(total))
            session.commit()

        orig_v_genes = set(
            re.findall('IGHV[^ ,]+', seq['record'][columns.v_gene])
        )
        orig_j_genes = set(
            re.findall('IGHJ[^ ,]+', seq['record'][columns.j_gene])
        )
        if remaps is not None:
            remapped_j_genes = set([])
            for j in orig_j_genes:
                for remap_from, remap_to in remaps.iteritems():
                    if j.startswith(remap_from):
                        remapped_j_genes.add(remap_to)
                        break
                else:
                    remapped_j_genes.add(j)
            orig_j_genes = remapped_j_genes

        v_genes = filter(lambda v: v in v_germlines, orig_v_genes)
        j_genes = filter(lambda j: j in j_germlines, orig_j_genes)

        vdj = VDJSequence(
            seq['seq_ids'], seq['record'][columns.full_sequence], v_germlines,
            j_germlines, force_vs=v_genes, force_js=j_genes
        )
        try:
            if len(v_genes) == 0:
                raise AlignmentException('No valid V germline for {}'.format(
                    ','.join(sorted(orig_v_genes))
                ))
            if len(j_genes) == 0:
                raise AlignmentException('No valid J germline for {}'.format(
                    ','.join(sorted(orig_j_genes))
                ))
            vdj.analyze()

            if vdj.sequence in aligned_seqs:
                aligned_seqs[vdj.sequence].ids += vdj.ids
            else:
                aligned_seqs[vdj.sequence] = vdj
        except AlignmentException as e:
            add_as_noresult(session, vdj, sample, str(e))
            missed += 1
    logger.info('Aligned {} / {} sequences'.format(total - missed + 1, total))

    logger.info('Collapsing ambiguous character sequences')
    if len(aligned_seqs) > 0:
        avg_mut = sum(
            [v.mutation_fraction for v in aligned_seqs.values()]
        ) / float(len(aligned_seqs))
        avg_len = sum(
            [v.v_length for v in aligned_seqs.values()]
        ) / float(len(aligned_seqs))
        sample.v_ties_mutations = avg_mut
        sample.v_ties_len = avg_len
        if columns.ties:
            add_uniques(session, sample, aligned_seqs.values(),
                        realign_mut=avg_mut, realign_len=avg_len,
                        trim_to=columns.trim_to,
                        max_padding=columns.max_padding)
        else:
            add_uniques(session, sample, aligned_seqs.values())
    session.commit()


def run_import(session, args, remaps=None):
    v_germlines = VGermlines(args.v_germlines)
    j_germlines = JGermlines(args.j_germlines, args.upstream_of_cdr3,
                             args.anchor_len, args.min_anchor_len)

    study, new = funcs.get_or_create(session, Study, name=args.study_name)

    if new:
        logger.info('Created new study "{}"'.format(study.name))
        session.commit()

    sample, new = funcs.get_or_create(session, Sample, name=args.sample_name,
                                      study=study)
    if new:
        sample.date = args.date
        logger.info('Created new sample "{}"'.format(sample.name))
        for key in ('subset', 'tissue', 'disease', 'lab', 'experimenter',
                    'ig_class', 'v_primer', 'j_primer'):
            setattr(sample, key, vars(args).get(key, None))
        subject, new = funcs.get_or_create(
            session, Subject, study=study,
            identifier=args.subject)
        sample.subject = subject
        session.commit()
    else:
        loger.error('Sample "{}" already exists'.format(args.sample_name))
        return

    with open(args.input_file) as fh:
        read_file(session, fh, sample, v_germlines, j_germlines, args, remaps)
