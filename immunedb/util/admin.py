import getpass
import json
import pymysql
import os
import random
import re
import shlex
import string
import subprocess

import immunedb.common.config as config
from immunedb.util.log import logger


def _yn_prompt(prompt):
    while True:
        print('{} [Y/N]'.format(prompt), end=' ')
        v = input().lower()
        if v in ('y', 'n'):
            return v == 'y'


def _connect(host, user, admin_pass=''):
    return pymysql.connect(host=host, user=user, password=admin_pass,
                           cursorclass=pymysql.cursors.DictCursor)


def _get_root_connection(host, user, admin_pass=None):
    if admin_pass is None:
        try:
            return _connect(host, user)
        except Exception:
            logger.info('Failed connection with empty root password')
            admin_pass = getpass.getpass('MySQL password for ({}):'.format(
                user))
    return _connect(host, user, admin_pass)


def _create_user_if_not_exists(conn, host, user, password):
    with conn.cursor() as cursor:
        cursor.execute('SELECT host, user, authentication_string from '
                       'mysql.user WHERE user=%s and host=%s', (user, host))
        existing = cursor.fetchone()
        if existing is None:
            cursor.execute('CREATE USER \'{}\'@\'%\' IDENTIFIED BY '
                           '\'{}\''.format(user, password))
            return None
        return existing['authentication_string']


def _get_user_pass(conn, host, user, existing_password):
    with conn.cursor() as cursor:
        while True:
            db_pass = getpass.getpass()
            cursor.execute('SELECT PASSWORD(%s) as password', db_pass)
            if cursor.fetchone()['password'] != existing_password:
                logger.error('Password does not match.')
            else:
                logger.info('Correct password')
                return db_pass


def create(main_parser, args):
    if re.search(r'[^A-Za-z0-9_-]', args.db_name) is not None:
        main_parser.error('Database name must only contain letters, numbers, '
                          'dashes and underscores.')

    try:
        conn = _get_root_connection(args.db_host, args.admin_user,
                                    args.admin_pass)

        db_user = args.db_user or args.db_name
        if args.db_pass:
            db_pass = args.db_pass
        else:
            db_pass = ''.join(
                random.choice(string.ascii_uppercase + string.ascii_lowercase +
                              string.digits) for _ in range(10))

        with conn.cursor() as cursor:
            logger.info('Creating user "{}"'.format(db_user))
            existing_password = _create_user_if_not_exists(conn, '%', db_user,
                                                           db_pass)
            if existing_password is not None:
                if not args.db_pass:
                    logger.warning(
                        'User {} already exists.  To generate the '
                        'configuration file, you must enter it\'s '
                        'password.'.format(db_user)
                    )
                    db_pass = _get_user_pass(conn, args.db_host, db_user,
                                             existing_password)
                else:
                    db_pass = args.db_pass

            logger.info('Creating database "{}"'.format(args.db_name))
            cursor.execute('CREATE DATABASE {}'.format(args.db_name))

            cursor.execute(
                'GRANT ALL PRIVILEGES ON {}.* TO \'{}\'@\'%\''.format(
                    args.db_name, db_user))

        config_path = os.path.join(args.config_dir, '{}.json'.format(
            args.db_name))
        logger.info('Creating config at {}'.format(config_path))
        with open(config_path, 'w+') as fh:
            json.dump({
                'host': args.db_host,
                'database': args.db_name,
                'username': db_user,
                'password': db_pass
            }, fh, sort_keys=True, indent=4, separators=(',', ': '))

        logger.info('Initializing tables')
        config.init_db(config_path)
        logger.info('Success!')
        return True
    except Exception as e:
        logger.error(e)
        return False


def delete(main_parser, args):
    try:
        with open(args.db_config) as fh:
            db_config = json.load(fh)
        conn = _get_root_connection(db_config['host'], args.admin_user,
                                    args.admin_pass)
        with conn.cursor() as cursor:
            logger.info('Deleting database {}'.format(db_config['database']))
            cursor.execute('DROP DATABASE `{}`'.format(db_config['database']))
            if args.delete_user:
                logger.info('Deleting user {}'.format(db_config['username']))
                cursor.execute('DROP USER `{}`'.format(db_config['username']))
        return True
    except Exception as e:
        logger.error(e)
        return False


def backup(main_parser, args):
    with open(args.db_config) as fh:
        db_config = json.load(fh)
    with open(args.backup_path, 'w+') as fh:
        cmd = shlex.split(
            'mysqldump -h {} -u {} -p{} {}'.format(db_config['host'],
                                                   db_config['username'],
                                                   db_config['password'],
                                                   db_config['database']))
        proc = subprocess.Popen(cmd, stdout=fh)
        proc.communicate()

    return True


def restore(main_parser, args):
    with open(args.db_config) as fh:
        db_config = json.load(fh)
    with open(args.backup_path, 'r') as fh:
        cmd = shlex.split(
            'mysql -h {} -u {} -p{} {}'.format(db_config['host'],
                                               db_config['username'],
                                               db_config['password'],
                                               db_config['database']))
        proc = subprocess.Popen(cmd, stdin=fh, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        if stderr:
            logger.warning(stderr)
    return True
