import os

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


# V2 of multiprocessing
class _EndOfDataSentinel(object):
    pass


class EndOfDataException(Exception):
    pass


class SizedQueue(object):
    def __init__(self):
        self.queue = mp.Queue()

    def put(self, val):
        self.queue.put(val)

    def get(self):
        res = self.queue.get()
        if res == _EndOfDataSentinel:
            raise EndOfDataException()
        return res

    def __iter__(self):
        while True:
            try:
                yield self.get()
            except EndOfDataException:
                break


def _process_wrapper(process_func, in_queue, out_queue, **kwargs):
    logger.info('Starting worker with PID {}'.format(os.getpid()))
    while True:
        try:
            for d in process_func(in_queue.get(), **kwargs):
                out_queue.put(d)
        except Queue.Empty:
            continue
        except EndOfDataException:
            break
        except Exception:
            logging.warning('Error in process_func:\n{}'.format(
                traceback.format_exc()))
    logger.info('Stopping worker with PID {}'.format(os.getpid()))


def _agg_wrapper(aggregate_func, aggregate_queue, return_val, **kwargs):
    ret = aggregate_func(aggregate_queue, **kwargs)
    if ret:
        return_val.update(ret)


def process_data(generate_func_or_iter,
                 process_func,
                 aggregate_func,
                 nproc,
                 generate_args={},
                 process_args={},
                 aggregate_args={}):
    process_queue = SizedQueue()
    aggregate_queue = SizedQueue()
    return_val = mp.Manager().dict()

    workers = []
    for i in range(nproc):
        w = mp.Process(
            target=_process_wrapper,
            args=(process_func, process_queue, aggregate_queue),
            kwargs=process_args
        )
        w.start()
        workers.append(w)

    agg_p = mp.Process(
        target=_agg_wrapper,
        args=(aggregate_func, aggregate_queue, return_val,),
        kwargs=aggregate_args
    )
    agg_p.start()

    try:
        input_p = None
        iter(generate_func_or_iter)
        for v in generate_func_or_iter:
            process_queue.put(v)
    except TypeError:
        input_p = mp.Process(
            target=generate_func_or_iter,
            args=(process_queue,),
            kwargs=generate_args
        )
        input_p.start()

    if input_p:
        input_p.join()

    for n in range(nproc):
        process_queue.put(_EndOfDataSentinel)

    for i, w in enumerate(workers):
        w.join()
        logger.info('{} workers finished'.format(i + 1))
    aggregate_queue.put(_EndOfDataSentinel)

    logger.info('Waiting on aggregation')
    agg_p.join()
    logger.info('Done aggregation')

    return return_val.copy()
