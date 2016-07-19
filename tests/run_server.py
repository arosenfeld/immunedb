from gevent import monkey
monkey.patch_all()
from argparse import Namespace

from immunedb.api.rest_service import run_rest_service
import immunedb.common.config as config

session = config.init_db('test_db.json', as_maker=True)
run_rest_service(session, Namespace(
    port=8891,
    debug=True,
    allow_shutdown=True,
    rollbar_token=None,
    rollbar_env=None
))
