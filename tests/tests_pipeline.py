import json
import os
import unittest

from sldb.aggregation.clone_stats import run_clone_stats
import sldb.common.config as config
from sldb.common.models import Clone, CloneStats, Sequence, SequenceCollapse
from sldb.identification.identify import run_identify
from sldb.identification.local_align import run_fix_sequences
from sldb.aggregation.clones import run_clones
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
        self.local_align()
        self.collapse()
        self.clones()
        self.clone_stats()

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
        self.assertEqual(self.session.query(Sequence).count(), 1092)

        self._regression(
            'tests/data/post_identification.json',
            self.session.query(Sequence),
            'seq_id',
            ('bucket_hash', 'ai', 'seq_id', 'v_gene', 'j_gene',
             'num_gaps', 'pad_length', 'v_match', 'v_length', 'j_match',
             'j_length', 'pre_cdr3_length', 'pre_cdr3_match',
             'post_cdr3_length', 'post_cdr3_match', 'copy_number',
             'cdr3_num_nts', 'cdr3_num_nts', 'cdr3_nt', 'cdr3_aa',
             'sequence', 'quality', 'germline')
        )

    def local_align(self):
        checks = self.session.query(Sequence).filter(
            Sequence.probable_indel_or_misalign == 1
        ).all()
        run_fix_sequences(
            self.session,
            NamespaceMimic(
                v_germlines='tests/imgt_human_v.fasta',
                j_germlines='tests/imgt_human_j.fasta',
                upstream_of_cdr3=31,
                anchor_len=18,
                min_anchor_len=12,
                max_deletions=3,
                max_insertions=3
            )
        )
        self.assertEqual(
            self.session.query(Sequence).filter(
                Sequence.probable_indel_or_misalign == 1
            ).count(), 0
        )

    def collapse(self):
        run_collapse(
            self.session,
            NamespaceMimic(
                subject_ids=None
            )
        )
        self._regression(
            'tests/data/post_collapse.json',
            self.session.query(SequenceCollapse),
            'seq_ai',
            ('collapse_to_subject_sample_id',
             'collapse_to_subject_seq_ai',
             'collapse_to_subject_seq_id',
             'instances_in_subject',
             'copy_number_in_subject')
        )

    def clones(self):
        run_clones(
            self.session,
            NamespaceMimic(
                similarity=85,
                subject_ids=None,
                include_indels=False,
                exclude_partials=False,
                min_identity=0,
                min_copy=2,
                regen=False
            )
        )

        self._regression(
            'tests/data/post_clones_clones.json',
            self.session.query(Clone),
            'id',
            ('id', 'functional', 'v_gene', 'j_gene', 'insertions',
             'cdr3_nt', 'cdr3_num_nts', 'cdr3_aa', 'germline'),
        )
        self._regression(
            'tests/data/post_clones_assignment.json',
            self.session.query(Sequence),
            'seq_id',
            ('seq_id', 'clone_id'),
        )

    def clone_stats(self):
        run_clone_stats(
            self.session,
            NamespaceMimic(
                clone_ids=None,
                subject_ids=None,
                regen=False
            )
        )

        self._regression(
            'tests/data/post_clone_stats.json',
            self.session.query(CloneStats),
            'id',
            ('clone_id', 'sample_id', 'unique_cnt', 'total_cnt')
        )

    def _err(self, path, key, key_val, field, ref, val):
        return '{} @ {}={}: {} is {} instead of {}'.format(
            path, key, key_val, field, val, ref
        )

    def _generate(self, path, query, key, fields):
        print 'Generating regression for {}'.format(path)
        data = {}
        for record in query:
            data[getattr(record, key)] = {
                f: getattr(record, f) for f in fields
            }
        with open(path, 'w+') as fh:
            json.dump(data, fh, sort_keys=True, indent=4,
                      separators=(',', ': '))

    def _regression(self, path, query, key, fields):
        if os.getenv('GENERATE'):
            self._generate(path, query, key, fields)
        print 'Regression testing {}'.format(path)
        with open(path) as fh:
            checks = json.load(fh)
        for record in query:
            for fld, value in checks[str(getattr(record, key))].iteritems():
                self.assertEqual(
                    getattr(record, fld),
                    value,
                    self._err(path, key, getattr(record, key), fld, value,
                              getattr(record, fld))
                )

class NamespaceMimic(object):
    def __init__(self, **kwargs):
        self.nproc = 1
        self.db_config = CONFIG_PATH
        for k, v in kwargs.iteritems():
            setattr(self, k, v)
