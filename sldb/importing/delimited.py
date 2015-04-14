import csv

from sldb.common.models import Sample, Sequence, Study, Subject
from sldb.identification.v_genes import VGermlines
import sldb.util.funcs as funcs

IMPORT_HEADERS = {
    'study_name': 'The name of the study [Required]',

    'sample_name': 'The name of the sample [Required]',
    'subject': 'The name of the subject [Required]',

    'subset': 'The cell subset of the sample',
    'tissue': 'The tissue from which the sample was gathered',
    'disease': 'The disease present in the subject',
    'lab': 'The lab which gathered the sample',
    'experimenter': 'The individual who gathered the sample',

    'seq_id': 'A unique sequence identifier [Required]',
    'alignment': 'Read type for the sequence (R1, R2, or R1+R2) [Required]',
    'indel': 'A boolean indicating if the sequence has an indel',
    'v_gene': 'V-gene name [Required]',
    'j_gene': 'J-gene name [Required]',

    'in_frame': 'A boolean indicating if the sequence is in-frame [Required]',
    'functional': 'A boolean indicating if the sequence is functional '
                  '[Required]',
    'stop': 'A boolean indicating if the sequences contains any stop codons '
            '[Required]',
    'copy_number': 'The number of times the sequence occurred in the sample '
                   '[Required]',

    'sequence': 'The full, IMGT aligned sequence [Required]',
}


class ImportException(Exception):
    pass


class DelimitedImporter(object):
    def __init__(self, session, mappings, defaults):
        self._session = session
        self._mappings = mappings
        self._defaults = defaults
        self._cached_studies = {}
        self._cached_subjects = {}

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
                        'Header {} for field {} not found.'.format(header,
                            field_name))
                return None
            return self._defaults[header]
        return row[header]


    def process_file(self, fh, delimiter):
        for row in csv.DictReader(fh, delimiter=delimiter):
            study_name = self._get_value('study_name', row)
            sample_name = self._get_value('sample_name', row)

            cache_key = (study_name, sample_name)
            if cache_key in self._cached_studies:
                study, sample = self._cached_studies[cache_key]
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

                    subject_name = self._get_value('subject', row)
                    cache_key = (study_name, subject_name)
                    if cache_key in self._cached_subjects:
                        subject = self._cached_subjects[cache_key]
                    else:
                        subject, new = funcs.get_or_create(
                            self._session,
                            Subject,
                            study_id=study.id,
                            identifier=subject_name)
                        if new:
                            print 'Created new subject {}'.format(subject_name)
                    sample.subject = subject
                    self._session.flush()

                self._cached_studies[(study_name, sample_name)] = (study,
                        sample)
                # TODO: If not new, verify the rest of the fields are the same
        self._session.commit()


def _parse_map(lst):
    if lst is None:
        return {}
    mapping = {}
    for arg in lst:
        field, name = arg.split('=', 1)
        mapping[field] = name
    return mapping


def run_delimited_import(session, args):
    germlines = VGermlines(args.v_germlines, include_prepadded=True)
    mappings = _parse_map(args.mappings)
    defaults = _parse_map(args.defaults);

    importer = DelimitedImporter(session, mappings, defaults)

    for fn in args.files:
        with open(fn) as fh:
            importer.process_file(fh, args.delimiter)
