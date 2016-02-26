import json
import os
import unittest

import sldb.common.config as config
from sldb.common.models import (Clone, CloneStats, DuplicateSequence, NoResult,
                                Sample, SampleStats, Sequence,
                                SequenceCollapse)


DB_NAME = 'test_db'
CONFIG_PATH = '{}.json'.format(DB_NAME)


class RegressionTest(unittest.TestCase):
    def __init__(self, name, *args, **kwargs):
        super(RegressionTest, self).__init__(*args, **kwargs)
        self.name = name

    def setUp(self):
        self.session = config.init_db(CONFIG_PATH)

    def tearDown(self):
        self.session.close()

    def get_path(self, fn):
        return os.path.join('tests', 'data', 'regression', self.name, fn)

    def initial_regression(self):
        self._regression(
            self.get_path('post_identification.json'),
            self.session.query(Sequence),
            'seq_id',
            ('bucket_hash', 'ai', 'seq_id', 'v_gene', 'j_gene',
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
            ('name', 'subject_id', 'subset', 'tissue', 'ig_class', 'disease',
             'lab', 'experimenter', 'v_primer', 'j_primer',
             'v_ties_mutations', 'v_ties_len')
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
                agg_keys - checks_keys)
            self.assertEqual(checks.keys(), set(agg_keys))


class NamespaceMimic(object):
    def __init__(self, **kwargs):
        self.nproc = 1
        self.db_config = CONFIG_PATH
        for k, v in kwargs.iteritems():
            setattr(self, k, v)
