import multiprocessing as mp
import traceback
import Queue


class Worker(object):
    def _print(self, msg):
        print 'Worker {}: {}'.format(self._worker_id, msg)

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
            except Exception as ex:
                worker._print(
                    '[TASK ERROR] The task was not completed '
                    'because:\n{}'.format(traceback.format_exc())
                )
                self._task_queue.task_done()
        worker.cleanup()

    def num_tasks(self):
        return self._num_tasks


class TaskTerminator(object):
    pass


def _task_wrapper(func, read_queue, args=None):
    while True:
        val = read_queue.get()
        if val == TaskTerminator:
            return
        if args is not None:
            func(val, *args)
        else:
            func(val)

def setup_tasking(task_gen_func, worker_func, db_func, num_workers=None):
    if num_workers is None:
        num_workers = mp.cpu_count()
    queue_size = 2 * num_workers
    # Make the task and result queues
    task_queue = mp.Queue()
    result_queue = mp.Queue()

    # The task_creator generates tasks to process
    task_creator = mp.Process(target=task_gen_func, args=(task_queue,))
    task_creator.start()

    db_processor = mp.Process(target=_task_wrapper,
                              args=(db_func, result_queue))
    db_processor.start()


    # Make workers to process the tasks
    workers = []
    for i in range(0, num_workers):
        worker = mp.Process(target=_task_wrapper,
                            args=(worker_func, task_queue, (result_queue,)))
        worker.start()
        workers.append(worker)

    # Keep processing until no more tasks are added
    task_creator.join()
    # Add kill sentinels to the task queue
    for i in range(0, num_workers):
        task_queue.put(TaskTerminator)

    # Wait for all workers to complete
    for worker in workers:
        worker.join()
    result_queue.put(TaskTerminator)


if __name__ == '__main__':
    def task_gen_func(task_queue):
        for i in range(0, 100):
            task_queue.put(i)

    def worker_func(task, result_queue):
        result_queue.put((task, sum(range(0, task + 1))))

    def db_func(result):
        print 'Sum of 0..{} is {}'.format(*result)


    setup_tasking(task_gen_func, worker_func, db_func)
