import csv

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


def run_delimited_import(session, args):
    pass
    #germlines = VGermlines(args.v_germlines, include_prepadded=True)
