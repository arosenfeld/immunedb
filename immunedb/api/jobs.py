import os
import uuid
from functools import partial
import multiprocessing as mp

import logging
from immunedb.util.log import logger


class JobQueue:
    def __init__(self, temp_dir='/tmp/immunedb-jobs'):
        self.temp_dir = temp_dir
        if not os.path.exists(self.temp_dir):
            os.makedirs(self.temp_dir)
        self.files = []

    def get_path(self, uid, ext=''):
        return os.path.join(self.temp_dir, uid + ext)

    def _job_wrap(self, func, uid, **kwargs):
        logger.handlers = [
            logging.FileHandler(self.get_path(uid, '.log'))
        ]
        lock_fh = open(self.get_path(uid, '.LOCK'), 'wb+')
        with open(self.get_path(uid, '.zip'), 'wb+') as fh:
            val = func(**kwargs)
            fh.write(val)
        os.remove(lock_fh.name)

    def start_job(self, func, **kwargs):
        uid = str(uuid.uuid4())
        logger.info('Starting job with UUID {}'.format(uid))

        job_func = partial(self._job_wrap, func, uid)

        mp.Process(
            target=job_func,
            kwargs=kwargs
        ).start()
        self.files.extend([
            self.get_path(uid, '.log'),
            self.get_path(uid, '.zip'),
        ])

        return uid

    def job_complete(self, uid):
        if os.path.exists(self.get_path(uid, '.LOCK')):
            return False
        return True

    def get_log(self, uid):
        try:
            return open(self.get_path(uid, '.log'), 'r', buffering=1).read()
        except FileNotFoundError:
            return None

    def cleanup(self):
        for fn in self.files:
            try:
                os.remove(fn)
            except OSError as e:
                logger.warning('Could not remove {}: {}'.format(fn, e))
