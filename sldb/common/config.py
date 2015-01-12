import argparse
import json
import pkg_resources
import sys

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.schema import MetaData

from Bio import SeqIO

from sldb.common.settings import DATABASE_SETTINGS


class ConfigurationException(Exception):
    pass


def _create_engine(data):
    con_str = ('mysql+pymysql://{}:{}@{}/{}'
               '?charset=utf8&use_unicode=0').format(
        data['username'], data['password'],
        data['host'], data['database'])
    engine = create_engine(con_str, pool_recycle=3600)

    return engine, data['database']


def _get_germlines(path):
    if path.startswith('!'):
        try:
            path = pkg_resources.resource_filename(
                'sldb', 'data/germlines/imgt_human.fasta')
        except:
            raise ConfigurationException('Unknown built-in germline {}'.format(
                path[1:]))

    germs = {}
    try: 
        for record in SeqIO.parse(path, 'fasta'):
            germs[record.id] = str(record.seq)
    except:
        raise ConfigurationException('Unable to read V germlines at {}'.format(
            path))
    return germs


def get_base_arg_parser(desc):
    """Gets a base argument parser which requires a master and data
    configuration.

    :param str desc: The description provided by the argument parser

    :return: The argument parser object
    :rtype: ArgumentParser

    """
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('master_db_config', help='Path to master database'
                        'config')
    parser.add_argument('data_db_config', help='Path to data database config')

    return parser


def init_db(master_db_config, data_db_config, as_maker=False):
    """Initializes a session with the specified database database.

    :param str master_db_config: Path to master database config file
    :param str data_db_config: Path to data database config file
    :param bool as_maker: If ``True``, the returned object will be a session
        maker rather than an session

    :returns: A ``session`` or, if ``as_maker`` is set, a ``session_maker``

    """
    with open(master_db_config) as fh:
        master_config = json.load(fh)
    with open(data_db_config) as fh:
        data_config = json.load(fh)

    master_engine, master_name = _create_engine(master_config)
    data_engine, data_name = _create_engine(data_config)

    DATABASE_SETTINGS['master_metadata'] = MetaData(schema=master_name)
    DATABASE_SETTINGS['data_metadata'] = MetaData(schema=data_name)

    if 'v_germlines' not in data_config:
        raise ConfigurationException('No V germlines specified')
    DATABASE_SETTINGS['v_germlines'] = _get_germlines(
        data_config.get('v_germlines'))

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
