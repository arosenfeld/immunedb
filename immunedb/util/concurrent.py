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
class DataQueue(object):
    def __init__(self, num_writers):
        self.num_writers = num_writers
        self.queue = mp.Queue()
        self.num_writers_finished = mp.Value('i', 0, lock=True)
        self.size = mp.Value('i', 0, lock=True)

    def writer_finished(self):
        with self.num_writers_finished.get_lock():
            self.num_writers_finished.value += 1

    def put(self, data):
        with self.size.get_lock():
            self.size.value += 1
        self.queue.put(data)

    def put_all(self, data):
        for d in data:
            self.put(d)

    def get(self, log_progress=False):
        if (not self.size.value and
                self.num_writers_finished.value == self.num_writers):
            raise Queue.Empty
        with self.size.get_lock():
            self.size.value -= 1
            if (log_progress and self.size.value > 0 and
                    self.size.value % 1000 == 0):
                logger.info('{} remaining'.format(self.size.value))
        return self.queue.get()

    def __len__(self):
        return self.size.value


def _process_wrapper(process_func, in_queue, out_queue, log_progress,
                     **kwargs):
    i = 0
    while True:
        try:
            data = in_queue.get(log_progress=log_progress)
            i += 1
            out_queue.put(process_func(data, **kwargs))
        except Queue.Empty:
            break
        except Exception:
            logging.warning('Error in process_func:\n{}'.format(
                traceback.format_exc()))
    out_queue.writer_finished()


def process_data(generate_func_or_iter,
                 process_func,
                 aggregate_func,
                 nproc,
                 generate_args={},
                 process_args={},
                 aggregate_args={},
                 log_progress=False):
    process_queue = DataQueue(1)
    aggregate_queue = DataQueue(nproc)
    return_value = mp.Manager().dict()

    try:
        iter(generate_func_or_iter)
        process_queue.put_all(generate_func_or_iter)
        process_queue.writer_finished()
    except TypeError:
        mp.Process(
            target=generate_func_or_iter,
            args=(process_queue,),
            kwargs=generate_args
        ).start()

    for i in range(nproc):
        mp.Process(
            target=_process_wrapper,
            args=(process_func, process_queue, aggregate_queue, log_progress),
            kwargs=process_args
        ).start()
    aggregate_proc = mp.Process(
        target=aggregate_func,
        args=(aggregate_queue, return_value),
        kwargs=aggregate_args
    )
    aggregate_proc.start()
    aggregate_proc.join()
    return return_value
