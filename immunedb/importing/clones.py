import csv

from sqlalchemy.orm import joinedload

from immunedb.common.models import Clone, CloneStats, SampleStats, Sequence
from immunedb.aggregation.clones import generate_consensus, push_clone_ids
from immunedb.importing import ImportException


def generate_template(session, out_file):
    with open(out_file, 'w+') as fh:
        writer = csv.DictWriter(fh, delimiter='\t', fieldnames=[
            'ai',
            'seq_id',
            'sample_id',
            'subject_id',
            'subject',
            'sequence',
            'v_gene',
            'j_gene',
            'cdr3_aa',
            'cdr3_nt',
            'cdr3_num_nts',
            'copy_number',
            'copy_number_in_subject',
            'clone_id'
        ])

        writer.writeheader()
        template_seqs = session.query(Sequence).options(
            joinedload(Sequence.collapse)).order_by(Sequence.ai)
        for seq in template_seqs:
            if seq.collapse.copy_number_in_subject == 0:
                continue
            writer.writerow({
                'ai': seq.ai,
                'seq_id': seq.seq_id,
                'sample_id': seq.sample_id,
                'subject_id': seq.subject_id,
                'subject': seq.subject.identifier,
                'sequence': seq.sequence,
                'v_gene': seq.v_gene,
                'j_gene': seq.j_gene,
                'cdr3_aa': seq.cdr3_aa,
                'cdr3_nt': seq.cdr3_nt,
                'cdr3_num_nts': seq.cdr3_num_nts,
                'copy_number': seq.copy_number,
                'copy_number_in_subject': seq.collapse.copy_number_in_subject,
                'clone_id': 0
            })


def import_template(session, in_file, regen):
    if regen:
        session.query(Clone).delete()
        session.query(CloneStats).delete()
        session.query(SampleStats).delete()
        session.commit()

    seen_clones = {}
    with open(in_file) as fh:
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
                               v_gene=line['v_gene'],
                               j_gene=line['j_gene'],
                               cdr3_num_nts=int(line['cdr3_num_nts'])),
                'seqs': []
            })
            clone_inst = clone_info['clone']
            if (clone_inst.v_gene != line['v_gene'] or
                    clone_inst.j_gene != line['j_gene'] or
                    clone_inst.cdr3_num_nts != int(line['cdr3_num_nts']) or
                    clone_inst.subject_id != int(line['subject_id'])):
                raise ImportException(
                    'Sequence with ai {} was assigned to clone {} with '
                    'mismatching v_gene, j_gene, cdr3_num_nts, or '
                    'subject.'.format(line['ai'], line['clone_id']))
            clone_info['seqs'].append({
                'ai': int(line['ai']),
                'sample_id': int(line['sample_id'])
            })
        session.commit()

    db_clone_ids = set([])
    for cid, clone_info in sorted(seen_clones.items()):
        clone_inst = clone_info['clone']
        session.add(clone_inst)
        session.flush()
        db_clone_ids.add(clone_inst.id)

        to_update = [{
            'sample_id': s['sample_id'],
            'ai': s['ai'],
            'clone_id': clone_inst.id
        } for s in clone_info['seqs']]
        session.bulk_update_mappings(Sequence, to_update)
    session.commit()
    generate_consensus(session, db_clone_ids)
    push_clone_ids(session)
