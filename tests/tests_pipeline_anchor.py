from .regression import NamespaceMimic, BaseTest

from immunedb.identification.identify import run_identify


class TestPipelineAnchor(BaseTest.RegressionTest):
    def __init__(self, *args, **kwargs):
        super().__init__('pipeline_anchor', *args, **kwargs)

    def identification(self):
        run_identify(
            self.session,
            NamespaceMimic(
                v_germlines='tests/data/germlines/imgt_human_v.fasta',
                j_germlines='tests/data/germlines/imgt_human_j.fasta',
                upstream_of_cdr3=31,
                anchor_len=18,
                min_anchor_len=12,
                sample_dir='tests/data/identification',
                metadata=None,
                max_vties=50,
                min_similarity=0.60,
                trim=0,
                warn_existing=False,
                warn_missing=False,
                trim_to=None,
                max_padding=None,
                genotyping=False,
                ties=False,
            ),
        )
