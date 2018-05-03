import argparse
import json
import multiprocessing as mp

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from pymysql.cursors import SSCursor

from immunedb.common.models import Base


def get_base_arg_parser(desc='', multiproc=True, **kwargs):
    """Gets a base argument parser which requires database configuration.

    :param str desc: The description provided by the argument parser

    :return: The argument parser object
    :rtype: ArgumentParser

    """
    parser = argparse.ArgumentParser(
            description=desc,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            **kwargs)

    parser.add_argument('db_config', help='Path to database config')
    if multiproc:
        try:
            num_cpu = mp.cpu_count()
        except NotImplementedError:
            num_cpu = 4
        parser.add_argument('--nproc', default=num_cpu, type=int,
                            help='Number of subprocesses to run')

    return parser


def init_db(database_config, drop_all=False, as_maker=False, create=True):
    """Initializes a session with the specified database.

    :param str database_config: If ``from_dict`` is ``False``, the path to
        data database config file, otherwise a database config dictionary
    :param bool as_maker: If ``True``, the returned object will be a session
        maker rather than an session

    :returns: A ``session`` or, if ``as_maker`` is set, a ``session_maker``

    """

    if isinstance(database_config, str):
        with open(database_config) as fh:
            database_config = json.load(fh)

    conn = 'mysql+pymysql://{}:{}@{}/{}'.format(
        database_config['username'], database_config['password'],
        database_config['host'], database_config['database'])
    engine = create_engine(conn, pool_recycle=3600,
                           connect_args={'cursorclass': SSCursor})

    if drop_all:
        Base.metadata.drop_all(engine)

    if create:
        Base.metadata.create_all(engine)

    session = sessionmaker()
    session.configure(bind=engine)

    return session if as_maker else session()
