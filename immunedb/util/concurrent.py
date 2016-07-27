import multiprocessing as mp
import traceback
import Queue

import logging
from immunedb.util.log import logger


class Worker(object):
    def log(self, lvl, msg):
        logger.log(lvl, 'Worker {}: {}'.format(self._worker_id, msg))

    def info(self, msg):
        self.log(logging.INFO, msg)

    def error(self, msg):
        self.log(logging.ERROR, msg)

    def warning(self, msg):
        self.log(logging.WARNING, msg)

    def do_task(self, args):
        raise NotImplementedError

    def cleanup(self):
        pass


class TaskQueue(object):
    def __init__(self):
        self._task_queue = mp.JoinableQueue()
        self._num_tasks = 0
        self._workers = []

    def add_task(self, args):
        self._num_tasks += 1
        self._task_queue.put(args)

    def add_tasks(self, tasks):
        for task in tasks:
            self.add_task(task)

    def add_worker(self, worker):
        self._workers.append(
            mp.Process(
                target=self._func_wrap,
                args=(len(self._workers) + 1, worker)
            )
        )

    def start(self, block=True):
        for _ in self._workers:
            self.add_task(None)

        for worker in self._workers:
            worker.start()
        if block:
            self._task_queue.join()

    def _func_wrap(self, worker_id, worker):
        worker._worker_id = worker_id
        while True:
            try:
                try:
                    args = self._task_queue.get()
                except Queue.Empty:
                    break
                if args is None:
                    self._task_queue.task_done()
                    break
                else:
                    worker.do_task(args)
                    self._task_queue.task_done()
            except Exception:
                worker.error(
                    'The task was not completed because:\n{}'.format(
                        traceback.format_exc()))
                self._task_queue.task_done()
        worker.cleanup()

    def num_tasks(self):
        return self._num_tasks
