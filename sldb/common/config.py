import argparse
import json
import sys

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.schema import MetaData

from sldb.common.settings import DATABASE_SETTINGS


def _create_engine(config_path):
    with open(config_path) as fh:
        data = json.load(fh)

    con_str = ('mysql+pymysql://{}:{}@{}/{}'
               '?charset=utf8&use_unicode=0').format(
        data['username'], data['password'],
        data['host'], data['database'])
    engine = create_engine(con_str)

    return engine, data['database']


def get_base_arg_parser(desc):
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('master_db_config', help='Path to master database'
                        'config')
    parser.add_argument('data_db_config', help='Path to data database config')

    return parser


def init_db(master_db_config, data_db_config, as_maker=False):
    master_engine, master_name = _create_engine(master_db_config)
    data_engine, data_name = _create_engine(data_db_config)

    DATABASE_SETTINGS['master_metadata'] = MetaData(schema=master_name)
    DATABASE_SETTINGS['data_metadata'] = MetaData(schema=data_name)
    from sldb.common.models import *
    BaseMaster.metadata.create_all(master_engine)
    BaseData.metadata.create_all(data_engine)

    model_map = {
        Study: master_engine,
        Sample: master_engine,
        CloneGroup: master_engine,
        Subject: master_engine,

        SampleStats: data_engine,
        Sequence: data_engine,
        SequenceMapping: data_engine,
        DuplicateSequence: data_engine,
        Clone: data_engine,
        NoResult: data_engine,
    }

    session = sessionmaker(twophase=True)
    session.configure(binds=model_map)
    if not as_maker:
        session = session()

    return session
