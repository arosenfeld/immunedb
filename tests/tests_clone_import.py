import os

from .regression import BaseTest

import immunedb.common.config as config
from immunedb.common.models import Clone, Sequence
from immunedb.importing.clones import generate_template, import_template

from .regression import CONFIG_PATH


class TestCloneImport(BaseTest.BaseRegression):
    def __init__(self, *args, **kwargs):
        super(TestCloneImport, self).__init__('clone_import', *args, **kwargs)

    def setUp(self):
        self.session = config.init_db(CONFIG_PATH)

    def test_export(self):
        export_path = 'test_unassigned'
        generate_template(self.session, export_path)
        with open(export_path) as result:
            with open('tests/data/clone_import/unassigned.tsv') as expected:
                try:
                    assert result.read().strip() == expected.read().strip()
                finally:
                    os.remove(export_path)

    def test_import(self):
        import_template(
            self.session,
            'tests/data/clone_import/assigned.tsv',
            True
        )

        self.regression(
            self.get_path('post_clone_import_clones.json'),
            self.session.query(Clone),
            'id',
            ('id', 'functional', 'v_gene', 'j_gene', '_insertions',
                '_deletions', 'cdr3_nt', 'cdr3_num_nts', 'cdr3_aa',
                'germline'),
        )
        self.regression(
            self.get_path('post_clone_import_assignment.json'),
            self.session.query(Sequence),
            'seq_id',
            ('seq_id', 'clone_id'),
        )
