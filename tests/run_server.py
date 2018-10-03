from argparse import Namespace  # noqa: E402

from immunedb.api.rest_service import run_rest_service  # noqa: E402
import immunedb.common.config as config  # noqa: E402

run_rest_service(Namespace(
    db_config='test_db.json',
    nproc=1,
    port=8891,
    debug=True,
    allow_shutdown=True,
    rollbar_token=None,
    rollbar_env=None,
))
