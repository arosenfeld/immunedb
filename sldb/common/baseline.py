import csv
import os
import shlex
import subprocess

from sqlalchemy.sql import distinct

from sldb.common.models import Clone, Sequence
from sldb.identification.v_genes import VGene
import sldb.common.config as config

TEST_FOCUSED = 1
TEST_LOCAL = 2
SPECIES_HUMAN = 1
SPECIES_MOUSE = 2
SUB_UNIFORM = 0
SUB_SMITH = 1
MUT_UNIFORM = 0
MUT_SHAPIRO = 1
CONSTANT_BOUNDARIES = [1, 26, 38, 55, 65, 104]

SEQ_CLONAL = 1
FIX_INDELS = 1


def get_selection(session, clone_id, script_path, samples=None,
                  temp_dir='/tmp',
                  test_type=TEST_FOCUSED,
                  species=SPECIES_HUMAN,
                  sub_model=SUB_UNIFORM,
                  mut_model=MUT_UNIFORM):
    clone = session.query(Clone).filter(Clone.id == clone_id).first()
    last_region = CONSTANT_BOUNDARIES[-1] + clone.cdr3_num_nts // 3
    boundaries = '{}:{}'.format(':'.join(map(str, CONSTANT_BOUNDARIES)),
                                last_region)
    unique_id = '_{}_{}'.format(clone_id, '_'.join(map(str, samples)))
    input_path = os.path.join(temp_dir, 'clone{}.fasta'.format(unique_id))
    out_path = os.path.join(temp_dir, 'output{}'.format(unique_id))
    read_path = os.path.join(temp_dir, 'output{}{}.txt'.format(unique_id,
        clone.id))

    _make_input_file(session, input_path, clone, samples)
    cmd = 'Rscript {} {} {} {} {} {} {} {} {} {} {}'.format(
        script_path, test_type, species,
        sub_model, mut_model, SEQ_CLONAL,
        FIX_INDELS, boundaries, input_path, out_path,
        clone.id)
    proc = subprocess.Popen(shlex.split(cmd),
                            cwd=os.path.dirname(script_path),
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    proc.communicate()

    with open(read_path) as fh:
        output = _parse_output(session, clone, fh)

    os.unlink(input_path)
    os.unlink(read_path)
    os.unlink(os.path.join(
        temp_dir, 'output{}{}.RData'.format(unique_id, clone.id)))

    return output


def _make_input_file(session, input_path, clone, samples):
    with open(input_path, 'w+') as fh:
        fh.write('>>>CLONE\n')
        fh.write('>>germline\n')
        germline = ''.join([
            clone.group.germline[:VGene.CDR3_OFFSET],
            clone.cdr3_nt,
            clone.group.germline[VGene.CDR3_OFFSET + clone.cdr3_num_nts:]
        ])
        fh.write('{}\n'.format(germline))
        query = session.query(
            distinct(Sequence.sequence_replaced).label('seq')
        ).filter(
            Sequence.clone_id == clone.id
        )
        if samples is not None:
            query = query.filter(Sequence.sample_id.in_(samples))

        for i, seq in enumerate(query):
            fh.write('>{}\n{}\n'.format(i, seq.seq))


def _parse_output(session, clone, fh):
    reader = csv.DictReader(fh, delimiter='\t')
    for row in reader:
        if row['Type'] == 'Sequence':
            del row['Type']
            del row['ID']
            row = {k: v.strip() for k, v in row.iteritems()}
            row = {
                k: v.strip() if v == 'NA' else float(v.strip()) for k, v in
                row.iteritems()
            }
            return row
