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
        print(f'{prompt} [Y/N]', end=' ')
        v = input().lower()
        if v in ('y', 'n'):
            return v == 'y'


def _connect(host, user, admin_pass=''):
    return pymysql.connect(
        host=host,
        user=user,
        password=admin_pass,
        cursorclass=pymysql.cursors.DictCursor,
    )


def _get_root_connection(host, user, admin_pass=None):
    if admin_pass is None:
        try:
            return _connect(host, user)
        except Exception:
            logger.info('Failed connection with empty root password')
            admin_pass = getpass.getpass(
                'MySQL password for ({}):'.format(user)
            )
    return _connect(host, user, admin_pass)


def _create_user_if_not_exists(conn, host, user, password):
    with conn.cursor() as cursor:
        cursor.execute(
            'SELECT host, user, authentication_string from '
            'mysql.user WHERE user=%s and host=%s',
            (user, host),
        )
        existing = cursor.fetchone()
        if existing is None:
            cursor.execute(
                'CREATE USER \'{}\'@\'%\' IDENTIFIED BY '
                '\'{}\''.format(user, password)
            )
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
        main_parser.error(
            'Database name must only contain letters, numbers, '
            'dashes and underscores.'
        )

    try:
        conn = _get_root_connection(
            args.db_host, args.admin_user, args.admin_pass
        )

        db_user = args.db_user or args.db_name
        if args.db_pass:
            db_pass = args.db_pass
        else:
            db_pass = ''.join(
                random.choice(
                    string.ascii_uppercase
                    + string.ascii_lowercase
                    + string.digits
                )
                for _ in range(10)
            )

        with conn.cursor() as cursor:
            logger.info(f'Creating user "{db_user}"')
            existing_password = _create_user_if_not_exists(
                conn, '%', db_user, db_pass
            )
            if existing_password is not None:
                if not args.db_pass:
                    logger.warning(
                        'User {} already exists.  To generate the '
                        'configuration file, you must enter it\'s '
                        'password.'.format(db_user)
                    )
                    db_pass = _get_user_pass(
                        conn, args.db_host, db_user, existing_password
                    )
                else:
                    db_pass = args.db_pass

            logger.info(f'Creating database "{args.db_name}"')
            cursor.execute(f'CREATE DATABASE {args.db_name}')

            cursor.execute(
                'GRANT ALL PRIVILEGES ON {}.* TO \'{}\'@\'%\''.format(
                    args.db_name, db_user
                )
            )

        config_path = os.path.join(
            args.config_dir, '{}.json'.format(args.db_name)
        )
        logger.info(f'Creating config at {config_path}')
        with open(config_path, 'w+') as fh:
            json.dump(
                {
                    'host': args.db_host,
                    'database': args.db_name,
                    'username': db_user,
                    'password': db_pass,
                },
                fh,
                sort_keys=True,
                indent=4,
                separators=(',', ': '),
            )

        logger.info('Initializing tables')
        config.init_db(config_path)
        _ensure_sequences_c_call(conn, args.db_name)
        logger.info('Success!')
        return True
    except Exception as e:
        logger.error(e)
        return False

# --- Added for AIRR constant-region support (c_call) ---
def _ensure_sequences_c_call(conn, db_name):
    """Ensure sequences.c_call exists (nullable) + index. Safe to run multiple times."""
    with conn.cursor() as cursor:
        cursor.execute(
            "SELECT COUNT(*) AS n FROM information_schema.COLUMNS "
            "WHERE TABLE_SCHEMA=%s AND TABLE_NAME='sequences' AND COLUMN_NAME='c_call';",
            (db_name,),
        )
        row = cursor.fetchone()
        exists = int(row["n"]) if isinstance(row, dict) else int(row[0])

        if exists == 0:
            logger.info("Adding sequences.c_call column")
            cursor.execute(f"USE `{db_name}`;")
            cursor.execute("ALTER TABLE sequences ADD COLUMN c_call VARCHAR(64) NULL;")
            cursor.execute("CREATE INDEX ix_sequences_c_call ON sequences (c_call);")
    conn.commit()
# --- end ---

def delete(main_parser, args):
    try:
        with open(args.db_config) as fh:
            db_config = json.load(fh)
        conn = _get_root_connection(
            db_config['host'], args.admin_user, args.admin_pass
        )
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
            'mysqldump -h {} -u {} -p{} {}'.format(
                db_config['host'],
                db_config['username'],
                db_config['password'],
                db_config['database'],
            )
        )
        proc = subprocess.Popen(cmd, stdout=fh)
        proc.communicate()

    return True


def restore(main_parser, args):
    with open(args.db_config) as fh:
        db_config = json.load(fh)
    with open(args.backup_path) as fh:
        cmd = shlex.split(
            'mysql -h {} -u {} -p{} {}'.format(
                db_config['host'],
                db_config['username'],
                db_config['password'],
                db_config['database'],
            )
        )
        proc = subprocess.Popen(
            cmd, stdin=fh, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        stdout, stderr = proc.communicate()
        if stderr:
            logger.warning(stderr)
    return True
