import csv
import re

from Bio import SeqIO

from sldb.common.models import CDR3_OFFSET, Sample, Sequence, Study, Subject
from sldb.identification.v_genes import get_common_seq
import sldb.util.funcs as funcs
import sldb.util.lookups as lookups

IMPORT_HEADERS = {
    'study_name': 'The name of the study [Required]',

    'sample_name': 'The name of the sample [Required]',
    'subject': 'The name of the subject [Required]',

    'paired': 'If the reads are paired end or not [Required]',

    'subset': 'The cell subset of the sample',
    'tissue': 'The tissue from which the sample was gathered',
    'disease': 'The disease present in the subject',
    'lab': 'The lab which gathered the sample',
    'experimenter': 'The individual who gathered the sample',

    'seq_id': 'A unique sequence identifier [Required]',
    'indel': 'A boolean indicating if the sequence has an indel',

    'v_gene': 'V-gene name [Required]',
    'v_length': 'Length of V gene EXCLUDING leading padding [Required]',
    'v_match': 'Number of nucleotides matching V-gene germline [Required]',

    'j_gene': 'J-gene name [Required]',
    'j_length': 'Length of J gene [Required]',
    'j_match': 'Number of nucleotides matching the J-gene germline [Required]',

    'in_frame': 'A boolean indicating if the sequence is in-frame [Required]',
    'functional': 'A boolean indicating if the sequence is functional '
                  '[Required]',
    'stop': 'A boolean indicating if the sequences contains any stop codons '
            '[Required]',
    'copy_number': 'The number of times the sequence occurred in the sample '
                   '[Required]',

    'sequence': 'The full, IMGT aligned sequence [Required]',
    'quality': 'Phred quality per-position in Sanger format.  Gapped '
               'positions should be spaces (Sanger format)',
    'cdr3_num_nts': 'The CDR3 number of nucleotides [Required]',
    'cdr3_nts': 'The CDR3 nucleotides [Required]',
    'cdr3_aas': 'The CDR3 amino-acids.  If not specified, the `cdr3_nts` will '
                'be converted to an amino-acid string with unknowns replaced '
                'with Xs',
}


class ImportException(Exception):
    pass


def _is_true(v):
    return v.upper() in ('T', 'TRUE')


class DelimitedImporter(object):
    def __init__(self, session, mappings, defaults, v_germlines, v_addition,
                 j_germlines, j_offset, fail_action):
        self._session = session
        self._mappings = mappings
        self._defaults = defaults
        self._v_germlines = v_germlines
        self._v_addition = v_addition
        self._j_germlines = j_germlines
        self._j_offset = j_offset
        self._fail_action = fail_action

        self._cached_studies = {}
        self._cached_subjects = {}
        self._cached_seqs = {}

    def _get_header_name(self, field_name):
        if field_name in self._mappings:
            return self._mappings[field_name]
        return field_name

    def _get_value(self, field_name, row, throw=True):
        header = self._get_header_name(field_name)
        if header not in row:
            if header not in self._defaults:
                if throw:
                    raise ImportException(
                        'Header {} for field {} not found.'.format(
                            header, field_name))
                return None
            return self._defaults[header]
        return row[header]

    def _get_models(self, row):
        study_name = self._get_value('study_name', row)
        sample_name = self._get_value('sample_name', row)
        subject_name = self._get_value('subject', row)

        sample_cache_key = (study_name, sample_name)
        if sample_cache_key in self._cached_studies:
            study, sample = self._cached_studies[sample_cache_key]
        else:
            study, new = funcs.get_or_create(self._session, Study,
                                             name=study_name)
            if new:
                print 'Created new study {}'.format(study_name)
                self._session.flush()

            sample, new = funcs.get_or_create(self._session, Sample,
                                              study_id=study.id,
                                              name=sample_name)
            if new:
                print 'Created new sample {}'.format(sample_name)
                sample.date = self._get_value('date', row)
                for field in ('subset', 'tissue', 'disease', 'lab',
                              'experimenter'):
                    setattr(sample, field,
                            self._get_value(field, row, throw=False))

                subject_cache_key = (study_name, subject_name)
                if subject_cache_key in self._cached_subjects:
                    subject = self._cached_subjects[subject_cache_key]
                else:
                    subject, new = funcs.get_or_create(
                        self._session,
                        Subject,
                        study_id=study.id,
                        identifier=subject_name)
                    if new:
                        print 'Created new subject {}'.format(subject_name)
                    self._cached_subjects[subject_cache_key] = subject

                sample.subject = subject
                self._session.flush()

            self._cached_studies[sample_cache_key] = (study, sample)
            # TODO: If not new, verify the rest of the fields are the same
        return study, sample

    def _process_sequence(self, row, study, sample):
        seq = self._get_value('sequence', row).upper().replace('.', '-')
        v_region = seq[:CDR3_OFFSET]
        pad_length = re.match('[-N]*', v_region).end() or 0
        # Some input uses gaps instead of N's, so replace them with N's
        seq = 'N' * pad_length + seq[pad_length:]

        vs = [g if g.startswith('IGHV') else 'IGHV{}'.format(g)
              for g in self._get_value('v_gene', row).split('|')]
        js = [g if g.startswith('IGHJ') else 'IGHJ{}'.format(g)
              for g in self._get_value('j_gene', row).split('|')]

        v_germline = get_common_seq(
            [self._v_germlines[v] for v in vs])[:CDR3_OFFSET]
        j_germline = get_common_seq(
            [self._j_germlines[j][-self._j_offset:] for j in js])

        germline = ''.join([
            v_germline,
            '-' * len(self._get_value('cdr3_nts', row)),
            j_germline
        ])

        if len(seq) < len(germline):
            seq += 'N' * (len(germline) - len(seq))

        seq_cache_key = (sample.id, seq)
        if seq_cache_key in self._cached_seqs:
            existing = self._cached_seqs[seq_cache_key]
        else:
            # Check for duplicate sequence
            existing = self._session.query(Sequence).filter(
                Sequence.sequence == seq,
                Sequence.sample_id == sample.id).first()

        if existing is not None:
            existing.copy_number += int(self._get_value('copy_number', row))
            self._session.flush()
            return
        cdr3_len = int(self._get_value('cdr3_num_nts', row))
        if cdr3_len == 0:
            cdr3_nt = cdr3_aa = ''
        else:
            cdr3_nt = self._get_value('cdr3_nts', row).upper()
            cdr3_aa = (self._get_value('cdr3_aas', row, throw=False) or
                       lookups.aas_from_nts(
                           self._get_value('cdr3_nts', row)))
        new_seq = Sequence(
            sample=sample,

            seq_id=self._get_value('seq_id', row),

            paired=_is_true(self._get_value('paired', row)),
            partial=pad_length > 0,
            probable_indel_or_misalign=_is_true(self._get_value('indel', row)),
            v_gene=self._get_value('v_gene', row),
            j_gene=self._get_value('j_gene', row),

            num_gaps=v_region[pad_length:].count('-'),
            pad_length=pad_length,

            v_match=self._get_value('v_match', row),
            v_length=int(self._get_value('v_length', row)),

            j_match=self._get_value('j_match', row),
            j_length=self._get_value('j_length', row),

            pre_cdr3_length=int(self._get_value('v_length', row)) -
            self._v_addition,
            pre_cdr3_match=self._get_value('v_match', row),
            post_cdr3_length=len(j_germline),
            post_cdr3_match=self._get_value('j_match', row),

            in_frame=_is_true(self._get_value('in_frame', row)),
            functional=_is_true(self._get_value('functional', row)),
            stop=_is_true(self._get_value('stop', row)),
            copy_number=int(self._get_value('copy_number', row)),

            cdr3_num_nts=cdr3_len,
            cdr3_nt=cdr3_nt,
            cdr3_aa=cdr3_aa,

            gap_method='IMGT',

            sequence=seq,
            quality=self._get_value('quality', row, throw=False),

            germline=germline,
        )

        self._cached_seqs[seq] = new_seq
        self._session.add(new_seq)

    def process_file(self, fh, delimiter, samples_to_collapse):
        for i, row in enumerate(csv.DictReader(fh, delimiter=delimiter)):
            try:
                study, sample = self._get_models(row)
                self._process_sequence(row, study, sample)
                if i % 1000 == 0 and i > 0:
                    print 'Processed {} sequences'.format(i)
                samples_to_collapse.add(sample.id)
            except Exception as ex:
                if self._fail_action != 'pass':
                    print ('[WARNING] Unable to process row #{}: '
                           'exception={}, msg={}').format(
                        i, str(ex.__class__.__name__), ex.message)
                    if self._fail_action == 'fail':
                        raise ex

        self._session.commit()


def _parse_map(lst):
    if lst is None:
        return {}
    mapping = {}
    for arg in lst:
        field, name = arg.split('=', 1)
        mapping[field] = name
    return mapping


def _get_germlines(path):
    germs = {}
    with open(path) as fh:
        for record in SeqIO.parse(fh, 'fasta'):
            germs[record.id] = str(record.seq)
    return germs


def run_delimited_import(session, args):
    mappings = _parse_map(args.mappings)
    defaults = _parse_map(args.defaults)

    v_germlines = _get_germlines(args.v_germlines)
    j_germlines = _get_germlines(args.j_germlines)

    importer = DelimitedImporter(session, mappings, defaults, v_germlines,
                                 args.v_addition, j_germlines, args.j_offset,
                                 args.fail_action)

    samples_to_collapse = set([])
    for fn in args.files:
        with open(fn) as fh:
            importer.process_file(fh, args.delimiter, samples_to_collapse)
