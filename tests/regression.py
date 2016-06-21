import json
import os
import unittest

import immunedb.common.config as config
from immunedb.common.models import (Clone, CloneStats, DuplicateSequence,
                                    NoResult, Sample, SampleStats, Sequence,
                                    SequenceCollapse)
from immunedb.identification.local_align import run_fix_sequences
from immunedb.aggregation.clones import run_clones
from immunedb.aggregation.collapse import run_collapse
from immunedb.aggregation.clone_stats import run_clone_stats
from immunedb.aggregation.sample_stats import run_sample_stats


DB_NAME = 'test_db'
CONFIG_PATH = '{}.json'.format(DB_NAME)


class BaseTest(object):
    class RegressionTest(unittest.TestCase):
        def __init__(self, name, *args, **kwargs):
            super(BaseTest.RegressionTest, self).__init__(*args, **kwargs)
            self.name = name

        def setUp(self):
            self.session = config.init_db(CONFIG_PATH, drop_all=True)

        def tearDown(self):
            self.session.close()

        def get_path(self, fn):
            return os.path.join('tests', 'data', 'regression', self.name, fn)

        def testAll(self):
            self.identification()
            self.initial_regression()
            self.local_align()
            self.collapse()
            self.clones()
            self.clone_stats()
            self.sample_stats()

        def initial_regression(self):
            self._regression(
                self.get_path('post_identification.json'),
                self.session.query(Sequence),
                'seq_id',
                ('ai', 'seq_id', 'v_gene', 'j_gene',
                 'num_gaps', 'pad_length', 'v_match', 'v_length', 'j_match',
                 'j_length', 'pre_cdr3_length', 'pre_cdr3_match',
                 'post_cdr3_length', 'post_cdr3_match', 'copy_number',
                 'cdr3_num_nts', 'cdr3_num_nts', 'cdr3_nt', 'cdr3_aa',
                 'sequence', 'quality', 'germline')
            )

            self._regression(
                self.get_path('post_identification_samples.json'),
                self.session.query(Sample),
                'id',
                ('name', 'subject_id', 'subset', 'tissue', 'ig_class',
                 'disease', 'lab', 'experimenter', 'v_primer', 'j_primer',
                 'v_ties_mutations', 'v_ties_len')
            )

        def local_align(self):
            run_fix_sequences(
                self.session,
                NamespaceMimic(
                    v_germlines='tests/data/germlines/imgt_human_v.fasta',
                    j_germlines='tests/data/germlines/imgt_human_j.fasta',
                    align_path=os.getenv('LL_PATH'),
                    min_similarity=60,
                    upstream_of_cdr3=31,
                    max_deletions=3,
                    max_insertions=3,
                    max_padding=None
                )
            )
            self.session.commit()

            for seq in self.session.query(Sequence.v_gene):
                self.assertEqual(seq.v_gene.count('IGHV'), 1)

            self._regression(
                self.get_path('post_local_align_seqs.json'),
                self.session.query(Sequence),
                'seq_id',
                ('ai', 'seq_id', 'v_gene', 'j_gene', 'num_gaps',
                 'pad_length', 'v_match', 'v_length', 'j_match', 'j_length',
                 'pre_cdr3_length', 'pre_cdr3_match', 'post_cdr3_length',
                 'post_cdr3_match', 'copy_number', 'cdr3_num_nts',
                 'cdr3_num_nts', 'cdr3_nt', 'cdr3_aa', 'sequence', 'quality',
                 'germline', 'insertions', 'deletions')
            )

            self._regression(
                self.get_path('post_local_align_nores.json'),
                self.session.query(NoResult),
                'seq_id',
                ('sequence', 'quality')
            )

            self._regression(
                self.get_path('post_local_align_dups.json'),
                self.session.query(DuplicateSequence),
                'seq_id',
                ('duplicate_seq_ai',)
            )

        def collapse(self):
            run_collapse(
                self.session,
                NamespaceMimic(
                    subject_ids=None
                )
            )
            self.session.commit()

            self._regression(
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
            run_clones(
                self.session,
                NamespaceMimic(
                    similarity=85,
                    subject_ids=None,
                    include_indels=False,
                    exclude_partials=False,
                    min_identity=0,
                    min_copy=2,
                    max_padding=None,
                    regen=False
                )
            )
            self.session.commit()

            self._regression(
                self.get_path('post_clones_clones.json'),
                self.session.query(Clone),
                'id',
                ('id', 'functional', 'v_gene', 'j_gene', 'insertions',
                 'cdr3_nt', 'cdr3_num_nts', 'cdr3_aa', 'germline',
                 'overall_unique_cnt', 'overall_total_cnt'),
            )
            self._regression(
                self.get_path('post_clones_assignment.json'),
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

            self._regression(
                self.get_path('post_clone_stats.json'),
                self.session.query(CloneStats),
                'id',
                ('clone_id', 'sample_id', 'subject_id', 'unique_cnt',
                 'total_cnt')
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

            self._regression(
                self.get_path('post_sample_stats.json'),
                self.session.query(SampleStats),
                ('sample_id', 'filter_type', 'outliers', 'full_reads'),
                ('sequence_cnt', 'in_frame_cnt', 'stop_cnt', 'functional_cnt',
                 'no_result_cnt')
            )

        def _err(self, path, key, key_val, field, ref, val):
            return '{} @ {}={}: {} is {} instead of {}'.format(
                path, key, key_val, field, val, ref
            )

        def _get_key(self, obj, key):
            if type(key) == str:
                return getattr(obj, key)
            return '-'.join([str(getattr(obj, k)) for k in key])

        def _generate(self, path, query, key, fields):
            print 'Generating regression for {}'.format(path)
            data = {}
            for record in query:
                data[self._get_key(record, key)] = {
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
            agg_keys = set([])
            for record in query:
                agg_key = str(self._get_key(record, key))
                agg_keys.add(agg_key)
                self.assertIn(agg_key, checks, '{} is in result but not in '
                              'checks for {}'.format(agg_key, path))
                for fld, value in checks[agg_key].iteritems():
                    self.assertEqual(
                        getattr(record, fld),
                        value,
                        self._err(path, key, self._get_key(record, key), fld,
                                  value, getattr(record, fld))
                    )
            check_keys = set(checks.keys())
            if check_keys != agg_keys:
                print 'Keys differ:'
                print '\tIn checks but not result: {}'.format(
                    check_keys - agg_keys)
                print '\tIn result but not checks: {}'.format(
                    agg_keys - check_keys)
                self.assertEqual(check_keys, set(agg_keys))


class NamespaceMimic(object):
    def __init__(self, **kwargs):
        self.nproc = 1
        self.db_config = CONFIG_PATH
        for k, v in kwargs.iteritems():
            setattr(self, k, v)
