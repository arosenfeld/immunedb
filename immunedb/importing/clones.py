import csv

from immunedb.common.models import Sequence, SequenceCollapse, Subject
from immunedb.aggregates.clones import generate_consensus
from immunedb.importing import ImportException

def generate_template(session, out_file):
    with open(out_file, 'w+') as fh:
        writer = csv.DictWriter(fh, delimiter='\t', fieldnames=[
            'ai',
            'seq_id',
            'subject_id',
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
        template_seqs = session.query(
            Sequence).join(SequenceCollapse).join(Subject)
        for seq in template_seqs:
            writer.writerow({
                'ai': seq.ai,
                'seq_id': seq.seq_id,
                'subject_id': seq.subject.id,
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
            clone_info = seen_clones.setdefault(line['clone_id'], {
                'clone': Clone(subject_id=int(line['subject_id']),
                               v_gene: line['v_gene'],
                               j_gene: line['j_gene'],
                               cdr3_num_nts: len(line['cdr3_nt'])),
                'seqs': []
            })
            clone_inst = clone_info['clone']
            if (clone_inst['v_gene'] != line['v_gene'] or
                    clone_inst['j_gene'] != line['j_gene'] or
                    clone_inst['cdr3_num_nts'] != len(line['cdr3_nt']),
                    clone_inst['subject_id'] == int(line['subject_id'])):
                raise ImportException(
                    'Sequence with ai {} was assigned to clone {} with '
                    'mismatching v_gene, j_gene, cdr3_num_nts, or '
                    'subject.'.format(line['ai'], line['clone_id']))
            clone_info['seqs'].append({
                'ai': line['ai'],
                'cdr3_nt': line['cdr3_nt']
            })


def create_clones(session, clone_mapping):
    db_clone_ids = set([])
    for clone_id, clone_info in clone_mapping.iteritems():
        clone_inst = clone_info['clone']
        session.add(clone_inst)
        session.flush()
        db_clone_ids.add(clone_inst.id)
    generate_consensus(session, db_clone_ids)
