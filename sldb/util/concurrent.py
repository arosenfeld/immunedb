import multiprocessing as mp
import Queue

class Worker(object):
    def do_task(self, worker_id, args):
        raise NotImplementedError

    def cleanup(self, worker_id):
        raise NotImplementedError


class TaskQueue(object):
    def __init__(self):
        self._task_queue = mp.JoinableQueue()
        self._workers = []

    def add_task(self, args):
        self._task_queue.put(args)

    def add_worker(self, worker):
        self._workers.append(
            mp.Process(target=self._func_wrap, args=(len(self._workers),
                worker))
        )

    def start(self, block=True):
        for worker in self._workers:
            worker.start()
        if block:
            self._task_queue.join()

    def _func_wrap(self, worker_id, worker):
        while True:
            try:
                args = self._task_queue.get_nowait()
            except Queue.Empty:
                break
            worker.do_task(worker_id, args)
            self._task_queue.task_done()
        worker.cleanup(worker_id)
