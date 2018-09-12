from argparse import Namespace  # noqa: E402

from immunedb.api.rest_service import run_rest_service  # noqa: E402
import immunedb.common.config as config  # noqa: E402

session = config.init_db('test_db.json', as_maker=True)
run_rest_service(session, Namespace(
    port=8891,
    debug=True,
    allow_shutdown=True,
    rollbar_token=None,
    rollbar_env=None,
    server='wsgiref'
))
