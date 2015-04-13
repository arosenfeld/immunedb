import argparse
import json
import multiprocessing as mp
import pkg_resources
import sys

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy.schema import MetaData

from pymysql.cursors import SSCursor

from Bio import SeqIO

from sldb.common.settings import DATABASE_SETTINGS

allowed_read_types = ('R1', 'R2', 'R1+R2')


def _create_engine(data):
    con_str = ('mysql+pymysql://{}:{}@{}/{}'
               '?charset=utf8&use_unicode=0').format(
        data['username'], data['password'],
        data['host'], data['database'])
    engine = create_engine(
        con_str, pool_recycle=3600, connect_args={
            'cursorclass': SSCursor})

    return engine, data['database']


def get_base_arg_parser(desc, multiproc=True, **kwargs):
    """Gets a base argument parser which requires a master and data
    configuration.

    :param str desc: The description provided by the argument parser

    :return: The argument parser object
    :rtype: ArgumentParser

    """
    parser = argparse.ArgumentParser(description=desc, **kwargs)
    parser.add_argument('master_db_config', help='Path to master database'
                        'config')
    parser.add_argument('data_db_config', help='Path to data database config')
    if multiproc:
        try:
            num_cpu = mp.cpu_count()
        except NotImplementedError:
            num_cpu = 4
        parser.add_argument('--nproc', default=num_cpu, type=int, help='Number'
                            ' of subprocesses to run ' '(Default: {})'.format(
                                num_cpu))

    return parser


def init_db(master_db_config, data_db_config, as_maker=False, from_dict=False):
    """Initializes a session with the specified database database.

    :param str master_db_config: If ``from_dict`` is ``False``, the path to
        master database config file, otherwise a database config dictionary
    :param str data_db_config: If ``from_dict`` is ``False``, the path to
        data database config file, otherwise a database config dictionary
    :param bool as_maker: If ``True``, the returned object will be a session
        maker rather than an session

    :returns: A ``session`` or, if ``as_maker`` is set, a ``session_maker``

    """

    if not from_dict:
        with open(master_db_config) as fh:
            master_db_config = json.load(fh)
        with open(data_db_config) as fh:
            data_db_config = json.load(fh)

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
        CloneStats: data_engine,
        Sequence: data_engine,
        DuplicateSequence: data_engine,
        Clone: data_engine,
        NoResult: data_engine,
        ModificationLog: data_engine,
    }

    session = sessionmaker(twophase=True)
    session.configure(binds=model_map)
    if not as_maker:
        session = session()

    return session
