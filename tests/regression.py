import json
import os
import unittest

from sqlalchemy.sql import text

import immunedb.common.config as config
from immunedb.common.models import (Clone, CloneStats, NoResult, Sample,
                                    SampleMetadata, SampleStats,
                                    SelectionPressure, Sequence,
                                    SequenceCollapse)
from immunedb.identification.local_align import run_fix_sequences
from immunedb.aggregation.clones import run_clones
from immunedb.aggregation.collapse import run_collapse
from immunedb.aggregation.clone_stats import run_clone_stats
from immunedb.aggregation.sample_stats import run_sample_stats
from immunedb.common.baseline import run_selection_pressure
from immunedb.trees.clearcut import run_clearcut

from .trees import tree_compare


DB_NAME = 'test_db'
CONFIG_PATH = '{}.json'.format(DB_NAME)


class BaseTest(object):
    class BaseRegression(unittest.TestCase):
        def __init__(self, name, *args, **kwargs):
            super(BaseTest.BaseRegression, self).__init__(*args, **kwargs)
            self.name = name

        def err(self, path, key, key_val, field, ref, val):
            return '{} @ {}={}: {} is {} instead of {}'.format(
                path, key, key_val, field, val, ref
            )

        def get_path(self, fn):
            return os.path.join('tests', 'data', 'regression', self.name, fn)

        def get_key(self, obj, key):
            if type(key) == str:
                return getattr(obj, key)
            return '-'.join([str(getattr(obj, k)) for k in key])

        def generate(self, path, query, key, fields):
            print('Generating regression for {}'.format(path))
            data = {}
            for record in query:
                data[self.get_key(record, key)] = {
                    f: getattr(record, f) for f in fields
                }
            with open(path, 'w+') as fh:
                json.dump(data, fh, sort_keys=True, indent=4,
                          separators=(',', ': '))

        def regression(self, path, query, key, fields, comp=None):
            if os.getenv('GENERATE'):
                self.generate(path, query, key, fields)
            print('Regression testing {}'.format(path))
            with open(path) as fh:
                checks = json.load(fh)
            agg_keys = set([])
            for record in query:
                agg_key = str(self.get_key(record, key))
                agg_keys.add(agg_key)
                self.assertIn(agg_key, checks, '{} is in result but not in '
                              'checks for {}'.format(agg_key, path))
                for fld, value in checks[agg_key].items():
                    (comp or self.assertEqual)(
                        getattr(record, fld),
                        value,
                        self.err(path, key, self.get_key(record, key), fld,
                                 value, getattr(record, fld))
                    )
            check_keys = set(checks.keys())
            if check_keys != agg_keys:
                print('Keys differ:')
                print('\tIn checks but not result: {}'.format(
                    check_keys - agg_keys))
                print('\tIn result but not checks: {}'.format(
                    agg_keys - check_keys))
                self.assertEqual(check_keys, set(agg_keys))

    class RegressionTest(BaseRegression):
        def __init__(self, name, *args, **kwargs):
            super(BaseTest.RegressionTest, self).__init__(name, *args,
                                                          **kwargs)

        def setUp(self):
            self.session_maker = config.init_db(CONFIG_PATH, drop_all=True,
                                                as_maker=True)
            self.session = self.session_maker()

        def tearDown(self):
            self.session.close()

        def testAll(self):
            self.identification()
            self.initial_regression()
            self.local_align()
            self.collapse()
            self.clones()
            self.clone_stats()
            self.sample_stats()
            self.trees()
            self.selection()

        def initial_regression(self):
            self.regression(
                self.get_path('post_identification.json'),
                self.session.query(Sequence),
                'seq_id',
                ('ai', 'seq_id', 'v_gene', 'j_gene',
                 'num_gaps', 'seq_start', 'v_match', 'v_length', 'j_match',
                 'j_length', 'pre_cdr3_length', 'pre_cdr3_match',
                 'post_cdr3_length', 'post_cdr3_match', 'copy_number',
                 'cdr3_num_nts', 'cdr3_num_nts', 'cdr3_nt', 'cdr3_aa',
                 'sequence', 'quality', 'germline')
            )

            self.regression(
                self.get_path('post_identification_samples.json'),
                self.session.query(Sample),
                'id',
                ('name', 'v_ties_mutations', 'v_ties_len')
            )

            self.regression(
                self.get_path('post_identification_sample_metadata.json'),
                self.session.query(SampleMetadata),
                ('sample_id', 'key'),
                ('value',)
            )

        def local_align(self):
            run_fix_sequences(
                self.session,
                NamespaceMimic(
                    v_germlines='tests/data/germlines/imgt_human_v.fasta',
                    j_germlines='tests/data/germlines/imgt_human_j.fasta',
                    temp='/tmp',
                    upstream_of_cdr3=31,
                    max_deletions=5,
                    max_insertions=5,
                    sample_ids=None
                )
            )
            self.session.commit()

            for seq in self.session.query(Sequence.v_gene):
                self.assertEqual(seq.v_gene.count('IGHV'), 1)

            self.regression(
                self.get_path('post_local_align_seqs.json'),
                self.session.query(Sequence),
                'seq_id',
                ('ai', 'seq_id', 'v_gene', 'j_gene', 'num_gaps',
                 'seq_start', 'v_match', 'v_length', 'j_match', 'j_length',
                 'pre_cdr3_length', 'pre_cdr3_match', 'post_cdr3_length',
                 'post_cdr3_match', 'copy_number', 'cdr3_num_nts',
                 'cdr3_num_nts', 'cdr3_nt', 'cdr3_aa', 'sequence', 'quality',
                 'germline', '_insertions', '_deletions')
            )

            self.regression(
                self.get_path('post_local_align_nores.json'),
                self.session.query(NoResult),
                'seq_id',
                ('sequence', 'quality')
            )

        def collapse(self):
            run_collapse(
                self.session,
                NamespaceMimic(
                    subject_ids=None
                )
            )
            self.session.commit()

            self.regression(
                self.get_path('post_collapse.json'),
                self.session.query(SequenceCollapse),
                'seq_ai',
                ('collapse_to_subject_sample_id',
                 'collapse_to_subject_seq_ai',
                 'collapse_to_subject_seq_id',
                 'instances_in_subject',
                 'copy_number_in_subject')
            )

        def clones(self):
            for method in ('cluster', 'similarity'):
                self.session.connection(mapper=Clone).execute(text('''
                    DELETE FROM clones
                '''))
                self.session.connection(mapper=Clone).execute(text('''
                    ALTER TABLE clones AUTO_INCREMENT=1
                '''))
                self.session.commit()
                run_clones(
                    self.session,
                    NamespaceMimic(
                        method=method,
                        subject_ids=None,
                        level='aa',
                        similarity=.85,
                        exclude_partials=False,
                        min_copy=2,
                        max_padding=None,
                        skip_regen=False,
                        gene=None,
                        reduce_difference=4,
                        skip_subclones=False,
                    )
                )
                self.session.commit()

                self.regression(
                    self.get_path('post_clones_clones_{}.json'.format(method)),
                    self.session.query(Clone),
                    'id',
                    ('id', 'functional', 'v_gene', 'j_gene', '_insertions',
                     '_deletions', 'cdr3_nt', 'cdr3_num_nts', 'cdr3_aa',
                     'germline', 'parent_id'),
                )
                self.regression(
                    self.get_path('post_clones_assignment_{}.json'.format(
                        method)),
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
            self.session.commit()

            self.regression(
                self.get_path('post_clone_stats.json'),
                self.session.query(CloneStats),
                'id',
                ('clone_id', 'sample_id', 'subject_id', 'functional',
                    'unique_cnt', 'total_cnt', 'avg_v_identity',
                    'top_copy_seq_ai')
            )

            self.regression(
                self.get_path('post_clone_stats_clones.json'),
                self.session.query(Clone),
                'id',
                ('id', 'overall_unique_cnt', 'overall_total_cnt',
                 'overall_instance_cnt'),
            )

        def sample_stats(self):
            run_sample_stats(
                self.session,
                NamespaceMimic(
                    sample_ids=None,
                    force=False
                )
            )
            self.session.commit()

            self.regression(
                self.get_path('post_sample_stats.json'),
                self.session.query(SampleStats),
                ('sample_id', 'filter_type', 'outliers', 'full_reads'),
                ('sequence_cnt', 'in_frame_cnt', 'stop_cnt', 'functional_cnt',
                 'no_result_cnt')
            )

        def trees(self):
            if not os.environ.get('TEST_TREES'):
                print('TEST_TREES not set.  Not testing trees')
                return
            run_clearcut(
                self.session,
                NamespaceMimic(
                    force=False,
                    clone_ids=None,
                    subject_ids=None,
                    temp='/tmp',
                    min_count=1,
                    min_seq_copies=0,
                    min_samples=1,
                    exclude_stops=False
                )
            )
            self.session.commit()

            self.regression(
                self.get_path('trees.json'),
                self.session.query(Clone),
                'id',
                ('tree',),
                comp=tree_compare
            )

        def selection(self):
            baseline_path = os.environ.get('BASELINE_PATH')

            if not baseline_path:
                print('Baseline path not set.  Skipping selection pressure '
                      'tests.')
                return
            run_selection_pressure(
                self.session,
                NamespaceMimic(
                    baseline_path=baseline_path,
                    # NOTE: Limit to 5 clones for speed sake
                    clone_ids=range(1, 6),
                    subject_ids=None,
                    regen=False,
                    temp='/tmp',
                    thresholds=['1', '2', '85%'],
                )
            )
            self.session.commit()

            self.regression(
                self.get_path('selection_pressure.json'),
                self.session.query(SelectionPressure),
                ('clone_id', 'sample_id', 'threshold'),
                ('expected_fwr_s', 'expected_cdr_s', 'expected_fwr_r',
                 'expected_cdr_r', 'observed_fwr_s', 'observed_cdr_s',
                 'observed_fwr_r', 'observed_cdr_r', 'sigma_fwr', 'sigma_cdr',
                 'sigma_fwr_cilower', 'sigma_fwr_ciupper', 'sigma_cdr_cilower',
                 'sigma_cdr_ciupper', 'sigma_p_fwr', 'sigma_p_cdr')
            )


class NamespaceMimic(object):
    def __init__(self, **kwargs):
        self.nproc = 1
        self.db_config = CONFIG_PATH
        for k, v in kwargs.items():
            setattr(self, k, v)
