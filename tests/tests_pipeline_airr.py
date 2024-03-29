from .regression import BaseTest, NamespaceMimic

from immunedb.importing.alignments import import_alignments


class TestPipelineAIRR(BaseTest.RegressionTest):
    def __init__(self, *args, **kwargs):
        super().__init__('pipeline_airr', *args, **kwargs)

    def identification(self):
        import_alignments(
            self.session,
            NamespaceMimic(
                v_germlines='tests/data/germlines/imgt_human_v.fasta',
                j_germlines='tests/data/germlines/imgt_human_j.fasta',
                sample_dir='tests/data/identification_import',
                metadata=None,
                max_vties=50,
                min_similarity=0.60,
                warn_existing=False,
                warn_missing=False,
                trim_to=None,
                max_padding=None,
                format='airr',
            ),
        )
        self.session.commit()

    def local_align(self):
        pass
