import csv as csv
import re

from Bio.Seq import Seq

from sldb.common.models import Sample, Study, Subject, Sequence
from sldb.identification.v_genes import VGermlines
from sldb.identification.vdj_sequence import VDJSequence
import sldb.util.funcs as funcs


class ImportException(Exception):
    pass


def _setup_import(session, args):
    study, new = funcs.get_or_create(session, Study,
                                     name=args.study_name)
    if new:
        print 'Created new study "{}"'.format(study.name)

    subject, new = funcs.get_or_create(session, Subject,
                                       study=study,
                                       identifier=args.subject)
    if new:
        print 'Created new subject "{}"'.format(subject.identifier)

    sample, new = funcs.get_or_create(session, Sample,
                                      study=study,
                                      name=args.sample_name)
    if new:
        print 'Created new sample "{}"'.format(sample.name)
        sample.date = args.date
        sample.subject = subject
        sample.subset = args.subset
        sample.tissue = args.tissue
        sample.disease = args.disease
        sample.lab = args.lab
        sample.experimenter = args.experimenter
    elif sample.subject != subject:
            raise ImportException('Tried to use existing sample with '
                                  'same name and different subject')
    session.commit()

    return sample


def _process_file(session, sample, reader, germlines, id_col, read_type, v_col,
                  j_col, identity_col, seq_col):
    # Sums of statistics for retroactive v-tie calculation
    lengths_sum = 0
    mutations_sum = 0
    # All VDJs assigned for this sample, keyed by raw sequence
    vdjs = {}
    # All probable duplicates based on keys of vdjs dictionary
    dups = {}

    for i, record in enumerate(reader):
        seq_id = record[id_col]
        identity = float(record[identity_col])
        sequence = Seq(record[seq_col].upper())

        if record[v_col] is None or record[j_col] is None:
            session.add(NoResult(
                seq_id=seq_id, sample=sample, sequence=sequence)
        v_gene = map(lambda g: 'IGHV{}'.format(g[0]),
                     re.findall('(\d+-\d+(\*\d+)?)', record[v_col]))
        j_gene = map(lambda g: 'IGHJ{}'.format(g[0]),
                     re.findall('(\d+)(\*\d+)?', record[j_col]))
        v_gene = filter(lambda v: v in germlines, v_gene)
        if len(v_gene) == 0 or len(j_gene) == 0:
            session.add(NoResult(
                seq_id=seq_id, sample=sample, sequence=sequence)

        # Key the vdjs dictionary by the unmodified sequence
        key = str(record.seq)
        if key in vdjs:
            # If this exact sequence, without padding or gaps, has been
            # assigned a V and J, bump the copy number of that
            # VDJSequence instance and add this seq_id as a duplicate.
            vdj = vdjs[key]
            vdj.copy_number += 1
            if vdj.id not in dups:
                dups[vdj.id] = []
            dups[vdj.id].append(DuplicateSequence(
                duplicate_seq_id=vdjs[key].id,
                sample_id=sample.id,
                seq_id=record.description))
        else:
            # This is the first instance of this exact sequence, so align it
            # and identify it's V and J
            vdj = VDJSequence(seq_id, sequence, read_type == 'R1+R2',
                              germlines, force_vs=v_gene)


        if i >= 10000:
            break


def run_csv_import(session, args):
    germlines = VGermlines(args.v_germlines)
    with open(args.csv_file) as fh:
        try:
            sample = _setup_import(session, args)
        except ImportException as ex:
            print '[ERROR] Unable setup import: \n\t{}'.format(
                ex.message)
            return
        reader = csv.DictReader(fh, delimiter=args.delim)
        _process_file(session, sample, reader, germlines,
                      id_col=args.id_col,
                      read_type=args.read_type,
                      v_col=args.v_col,
                      j_col=args.j_col,
                      identity_col=args.identity_col,
                      seq_col=args.seq_col)
