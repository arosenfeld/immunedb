import functools
import multiprocessing as mp
import traceback
import logging
import time

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
        for i, task in enumerate(tasks):
            self.add_task(task)

    def add_worker(self, worker):
        self._workers.append(
            mp.Process(
                target=self._func_wrap,
                args=(len(self._workers) + 1, worker)
            )
        )

    def add_workers(self, count, worker_cls, *args, **kwargs):
        for _ in range(count):
            self.add_worker(worker_cls(*args, **kwargs))

    def join(self):
        self._task_queue.join()

    def signal_end(self):
        for _ in self._workers:
            self.add_task(None)

    def start(self, block=True):
        if block:
            self.signal_end()

        for worker in self._workers:
            worker.start()
        if block:
            self.join()

    def _func_wrap(self, worker_id, worker):
        worker._worker_id = worker_id
        while True:
            try:
                args = self._task_queue.get()
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


def subcaller(func, data, i):
    return func(data[i])


# V2 of multiprocessing
def process_data(input_data, process_func, aggregate_func, nproc,
                 generate_args={}, process_args={}, aggregate_args={}):
    if callable(input_data):
        start = time.time()
        input_data = input_data(**generate_args)
        logger.info('Generate time: {}'.format(time.time() - start))

    with mp.Manager() as manager:
        proxy_data = manager.list(input_data)
        pool = mp.Pool(processes=nproc)
        f = functools.partial(
            subcaller,
            functools.partial(process_func, **process_args),
            proxy_data
        )
        start = time.time()
        logger.info('Waiting on pool {}'.format(process_func.__name__))

        res = [r for r in pool.map(f, range(len(proxy_data))) if r is not None]
        pool.close()
    logger.info('Pool done: {}'.format(time.time() - start))

    start = time.time()
    logger.info('Waiting on aggregation {}'.format(aggregate_func.__name__))
    ret = aggregate_func(res, **aggregate_args)
    logger.info('Done aggregation: {}'.format(time.time() - start))

    return ret
