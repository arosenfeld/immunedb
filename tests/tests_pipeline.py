import csv
import unittest

import sldb.common.config as config
from sldb.common.models import Sequence, SequenceCollapse
from sldb.identification.identify import run_identify
from sldb.aggregation.collapse import run_collapse

DB_NAME = 'test_db'
CONFIG_PATH = '{}.json'.format(DB_NAME)

class TestPipeline(unittest.TestCase):
    def setUp(self):
        self.session = config.init_db(CONFIG_PATH)

    def tearDown(self):
        self.session.close()

    def testAll(self):
        self.identification()
        self.collapse()

    def identification(self):
        run_identify(
            self.session,
            NamespaceMimic(
                v_germlines='tests/imgt_human_v.fasta',
                j_germlines='tests/imgt_human_j.fasta',
                upstream_of_cdr3=31,
                anchor_len=18,
                min_anchor_len=12,
                sample_dirs=['tests/data'],
                metadata=None,
                max_vties=50,
                min_similarity=60,
                trim=0,
                warn_existing=False,
            )
        )
        self.session.commit()
        self.assertTrue(
            self.session.query(Sequence).count() == 1092,
            self._err(
                'ALL', 'count', 1092, self.session.query(Sequence).count()
            )
        )
        with open('tests/data/post_identification.csv') as fh:
            checks = csv.DictReader(fh, delimiter='\t')
            for record in checks:
                seq = self.session.query(Sequence).filter(
                    Sequence.seq_id == record['seq_id']
                ).one()
                self._seq_checks(record, seq)

    def collapse(self):
        run_collapse(
            self.session,
            NamespaceMimic(
                subject_ids=None
            )
        )
        with open('tests/data/post_collapse.csv') as fh:
            checks = csv.DictReader(fh, delimiter='\t')
            simple_fields = (
                'collapse_to_subject_sample_id',
                'collapse_to_subject_seq_ai',
                'collapse_to_subject_seq_id',
                'instances_in_subject',
                'copy_number_in_subject'
            )
            for record in checks:
                collapse = self.session.query(SequenceCollapse).filter(
                    SequenceCollapse.sample_id == int(record['sample_id']),
                    SequenceCollapse.seq_ai == int(record['seq_ai'])
                ).one()
                for field in simple_fields:
                    self.assertEqual(
                        record[field],
                        str(getattr(collapse, field))
                    )

    def _seq_checks(self, reference, seq):
        simple_fields = (
            'bucket_hash', 'seq_id', 'v_gene', 'j_gene',
            'num_gaps', 'pad_length', 'v_match', 'v_length', 'j_match',
            'j_length', 'pre_cdr3_length', 'pre_cdr3_match',
            'post_cdr3_length', 'post_cdr3_match', 'copy_number',
            'cdr3_num_nts', 'cdr3_num_nts', 'cdr3_nt', 'cdr3_aa',
            'sequence', 'quality', 'germline'
        )
        for field in simple_fields:
            self.assertEqual(
                str(reference[field]), str(getattr(seq, field)),
                self._err(
                   seq.seq_id,
                   field,
                   str(reference[field]), str(getattr(seq, field))
                )
            )

        bool_fields = (
            'paired', 'partial', 'probable_indel_or_misalign', 'in_frame',
            'functional', 'stop'
        )
        for field in bool_fields:
            self.assertEqual(reference[field] == '1', getattr(seq, field))

    def _err(self, seq_id, field, ref, val):
        return '{}: {} should equal {}; instead {}'.format(
            seq_id,
            field,
            ref,
            val
        )


class NamespaceMimic(object):
    def __init__(self, **kwargs):
        self.nproc = 1
        self.db_config = CONFIG_PATH
        for k, v in kwargs.iteritems():
            setattr(self, k, v)
