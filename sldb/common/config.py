import argparse
import json

from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine

from models import *


def _create_engine(config_path):
    with open(config_path) as fh:
        data = json.load(fh)

    con_str = ('mysql://{}:{}@{}/{}?charset=utf8&use_unicode=0').format(
        data['username'], data['password'],
        data['host'], data['database'])
    engine = create_engine(con_str)

    return engine


def _create_tables(engine, tables):
    Base.metadata.create_all(engine, tables=tables)
    Base.metadata.bind = engine
    return sessionmaker(bind=engine)()


def get_base_arg_parser(desc):
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('master_db_config', help='Path to master database'
                        'config')
    parser.add_argument('data_db_config', help='Path to data database config')

    return parser


def get_session(parser):
    args = parser.parse_args()
    master_engine = _create_engine(args.master_db_config)
    data_engine = _create_engine(args.data_db_config)

    BaseMaster.metadata.create_all(master_engine)
    BaseData.metadata.create_all(data_engine)

    session = sessionmaker(twophase=True)
    session.configure(binds={
        Study: master_engine,
        Sample: master_engine,
        Clone: master_engine,
        Subject: master_engine,

        SampleStats: data_engine,
        Sequence: data_engine,
        DuplicateSequence: data_engine,
        NoResult: data_engine,
        CloneFrequency: data_engine
    })
    session = session()

    return session

if __name__ == '__main__':
    get_session(get_base_arg_parser('TEST'))
