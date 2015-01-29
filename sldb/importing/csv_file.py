import csv as csv
import re

from Bio.Seq import Seq

from sldb.common.models import (DuplicateSequence, NoResult, Sample, Study,
                                Subject, Sequence)
import sldb.identification.identify as identify
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
    no_result = 0

    for i, record in enumerate(reader):
        if record[v_col] is None or record[j_col] is None:
            session.add(NoResult(
                seq_id=seq_id, sample=sample, sequence=''))
            no_result += 1
            continue
        sequence = Seq(record[seq_col].upper())

        seq_id = record[id_col]
        identity = float(record[identity_col])

        v_gene = map(lambda g: 'IGHV{}'.format(g[0]),
                     re.findall('(\d+-\d+(\*\d+)?)', record[v_col]))
        j_gene = map(lambda g: 'IGHJ{}'.format(g[0]),
                     re.findall('(\d+)(\*\d+)?', record[j_col]))
        v_gene = filter(lambda v: v in germlines, v_gene)
        if len(v_gene) == 0 or len(j_gene) == 0:
            session.add(NoResult(
                seq_id=seq_id, sample=sample, sequence=sequence))
            no_result += 1
            continue

        # Key the vdjs dictionary by the unmodified sequence
        key = str(sequence)
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
                seq_id=seq_id))
        else:
            # This is the first instance of this exact sequence
            vdj = VDJSequence(seq_id, sequence, read_type == 'R1+R2',
                              germlines, force_vs=v_gene)
            if vdj.v_gene is not None and vdj.j_gene is not None:
                # If the V and J are found, add it to the vdjs dictionary to
                # prevent future exact copies from being aligned
                lengths_sum += vdj.v_length
                mutations_sum += vdj.mutation_fraction
                vdjs[key] = vdj
                print 'SUCCESS'
            else:
                if vdj.v_gene is None:
                    print 'NO V'
                elif vdj.j_gene is None:
                    print 'NO J'
                # The V or J could not be found, so add it as a noresult
                session.add(NoResult(sample=sample,
                                     seq_id=seq_id,
                                     sequence=str(vdj.sequence)))
                no_result += 1
        if i > 10000:
            break

    '''
    if len(vdjs) == 0:
        print '\t\tNo sequences identified'
        return

    print '\tCalculating V-ties'
    avg_len = lengths_sum / float(len(vdjs))
    avg_mut = mutations_sum / float(len(vdjs))

    for i, (_, vdj) in enumerate(vdjs.iteritems()):
        if i > 0 and i % 1000 == 0:
            print '\t\tCommitted {}'.format(i)
            session.commit()
        # Align the sequence to a germline based on v_ties
        vdj.align_to_germline(avg_len, avg_mut)
        if vdj.v_gene is not None and vdj.j_gene is not None:
            # Add the sequence to the database
            identify.add_to_db(session, read_type, sample, vdj)
        else:
            # This is a rare condition, but some sequence after aligning to
            # V-ties the CDR3 becomes non-existent, and it is thrown out
            no_result += 1
            session.add(NoResult(sample=sample,
                                 seq_id=vdj.id,
                                 sequence=str(vdj.sequence)))
            if vdj.id in dups:
                # It's duplicate sequences must be added as noresults also
                for dup in dups[vdj.id]:
                    # Restore the original sequence by removing padding and
                    # gaps
                    session.add(NoResult(
                        sample=sample,
                        seq_id=dup.seq_id,
                        sequence=vdj.sequence.replace('-', '').strip('N')))

                del dups[vdj.id]
    session.commit()
    '''

    # Add the true duplicates to the database
    print '\tAdding duplicates'
    for i, dup_seqs in enumerate(dups.values()):
        if i > 0 and i % 1000 == 0:
            print '\t\tCommitted {}'.format(i)
            session.commit()
        session.add_all(dup_seqs)
    session.commit()


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
