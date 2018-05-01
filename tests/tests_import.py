from regression import BaseTest, NamespaceMimic

from immunedb.importing.delimited import run_import


class TestImport(BaseTest.RegressionTest):
    def __init__(self, *args, **kwargs):
        super(TestImport, self).__init__('import', *args, **kwargs)

    def identification(self):
        run_import(
            self.session,
            NamespaceMimic(
                v_germlines='tests/data/germlines/imgt_human_v.fasta',
                j_germlines='tests/data/germlines/imgt_human_j.fasta',
                upstream_of_cdr3=31,
                anchor_len=18,
                min_anchor_len=12,
                sample_dir='tests/data/identification_import',
                metadata=None,
                max_vties=50,
                min_similarity=.60,
                warn_existing=False,
                warn_missing=False,
                trim=0,
                trim_to=None,
                max_padding=None,
            )
        )
        self.session.commit()
