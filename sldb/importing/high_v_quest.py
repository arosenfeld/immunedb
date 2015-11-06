import csv
import re

from sldb.identification import (add_as_noresult, add_as_sequence, add_uniques,
                                 AlignmentException)
from sldb.identification.vdj_sequence import VDJSequence
from sldb.identification.j_genes import JGermlines
from sldb.identification.v_genes import VGermlines
from sldb.common.models import NoResult, Sample, Study, Subject
import sldb.util.funcs as funcs

class ImportException(Exception):
    pass

def _collapse_seqs(session, sample, reader):
    seqs = {}
    for record in reader:
        if record['V-D-J-REGION'] is None:
            NoResult(seq_id=record['Sequence ID'], sample_id=sample.id)
            continue

        record['V-D-J-REGION'] = record['V-D-J-REGION'].replace(
            '.', '').upper()
        if record['V-D-J-REGION'] not in seqs:
            seqs[record['V-D-J-REGION']] = {
                'record': record,
                'seq_ids': []
            }
        seqs[record['V-D-J-REGION']]['seq_ids'].append(record['Sequence ID'])
    return seqs.values()

def read_file(session, handle, sample, v_germlines, j_germlines,
              paired):
    seqs = _collapse_seqs(session, sample, csv.DictReader(handle,
                          delimiter='\t'))

    aligned_seqs = {}
    missed = 0
    total = 0
    for total, seq in enumerate(seqs):
        if total > 0 and total % 1000 == 0:
            print 'Finished {}'.format(total)
            session.commit()
        v_genes = set(
            re.findall('IGHV[^ ]+', seq['record']['V-GENE and allele'])
        )
        j_genes = set(
            re.findall('IGHJ[^ ]+', seq['record']['J-GENE and allele'])
        )
        v_genes = filter(lambda v: v in v_germlines, v_genes)
        j_genes = filter(lambda j: j in j_germlines, j_genes)
        vdj = VDJSequence(
            seq['seq_ids'], seq['record']['V-D-J-REGION'], v_germlines,
            j_germlines, force_vs=v_genes, force_js=j_genes
        )
        try:
            if len(v_genes) == 0 or len(j_genes) == 0:
                raise AlignmentException('No V or J gene in input')
            vdj.analyze()
            vdj.align_to_germline()
            if vdj.sequence in aligned_seqs:
                aligned_seqs[vdj.sequence].ids += vdj.ids
            else:
                aligned_seqs[vdj.sequence] = vdj
        except AlignmentException as e:
            add_as_noresult(session, vdj, sample)
            missed += 1
    print 'Aligned {} / {} sequences'.format(total - missed, total)

    print 'Collapsing ambiguous character sequences'
    add_uniques(session, sample, aligned_seqs.values(), paired)
    session.commit()


def run_high_v_quest_import(session, args):
    v_germlines = VGermlines(args.v_germlines)
    j_germlines = JGermlines(args.j_germlines, args.upstream_of_cdr3,
                             args.anchor_len, args.min_anchor_len)

    study, new = funcs.get_or_create(session, Study, name=args.study_name)

    if new:
        print 'Created new study "{}"'.format(study.name)
        session.commit()

    sample, new = funcs.get_or_create(session, Sample, name=args.sample_name,
                                      study=study)
    if new:
        sample.date = args.date
        print 'Created new sample "{}"'.format(sample.name)
        for key in ('subset', 'tissue', 'disease', 'lab', 'experimenter',
                    'ig_class', 'v_primer', 'j_primer'):
            setattr(sample, key, vars(args).get(key, None))
        subject, new = funcs.get_or_create(
            session, Subject, study=study,
            identifier=args.subject)
        sample.subject = subject
        session.commit()
    else:
        print 'Sample already exists'
        return

    with open(args.gapped_nt_file) as fh:
        read_file(session, fh, sample, v_germlines, j_germlines,
                  not args.unpaired)
