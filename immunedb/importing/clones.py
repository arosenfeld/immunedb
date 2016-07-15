import csv

from immunedb.common.models import Sequence, SequenceCollapse, Subject
from immunedb.importing import ImportException

def generate_template(session, out_file):
    with open(out_file, 'w+') as fh:
        writer = csv.DictWriter(fh, delimiter='\t', fieldnames=[
            'ai',
            'seq_id',
            'subject',
            'sequence',
            'v_gene',
            'j_gene',
            'cdr3_aa',
            'cdr3_nt',
            'copy_number',
            'copy_number_in_subject',
            'clone_id'
        ])

        writer.writeheader()
        for seq in session.query(Sequence).join(SequenceCollapse).join(Subject):
            writer.writerow({
                'ai': seq.ai,
                'seq_id': seq.seq_id,
                'subject': seq.subject.identifier,
                'sequence': seq.sequence,
                'v_gene': seq.v_gene,
                'j_gene': seq.j_gene,
                'cdr3_aa': seq.cdr3_aa,
                'cdr3_nt': seq.cdr3_nt,
                'copy_number': seq.copy_number,
                'copy_number_in_subject': seq.collapse.copy_number_in_subject,
                'clone_id': 0
            })


def import_template(session, in_file):
    with open(in_file) as fh:
        seen_clones = {}
        reader = csv.DictReader(fh, delimiter='\t')
        if ('ai' not in reader.fieldnames or
                'clone_id' not in reader.fieldnames):
            raise ImportException(
                'Input file must have "ai" and "clone_id" fields')
        for line in reader:
            if line['clone_id'] in ('None', '0', ''):
                continue
            clone = seen_clones.setdefault(line['clone_id'], {
                'v_gene': line['v_gene'],
                'j_gene': line['j_gene'],
                'cdr3_num_nts': len(line['cdr3_nt']),
                'ais': []
            })
            print clone
            if (clone['v_gene'] != line['v_gene'] or
                    clone['j_gene'] != line['j_gene'] or
                    clone['cdr3_num_nts'] != len(line['cdr3_nt'])):
                raise ImportException('Sequence with ai {} was assigned to '
                                      'clone {} with mismatching v_gene, '
                                      'j_gene, or cdr3_num_nts'.format(
                                          line['ai'], line['clone_id']
                                      ))
            clone['ais'].append(line['ai'])
