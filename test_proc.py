import time
import random

from multiprocessing import Process, Queue, current_process, freeze_support



class Worker(Process):
    def __init__(self, queue):
        super(Worker, self).__init__()
        self.queue= queue

    def run(self):
        print 'Worker started'
        # do some initialization here

        print 'Computing things!'
        for data in iter( self.queue.get, None ):
            # Use data
            time.sleep(5)
            print data


request_queue = Queue()
for i in range(4):
    Worker( request_queue ).start()
for data in ['job%g'%x for x in range(10)]:
    request_queue.put( data )
# Sentinel objects to allow clean shutdown: 1 per worker.
for i in range(4):
    request_queue.put( None ) 








##
## Function run by worker processes
##
#
#def worker(input, output):
#    for func, args in iter(input.get, 'STOP'):
#        result = calculate(func, args)
#        output.put(result)
#
##
## Function used to calculate result
##
#
#def calculate(func, args):
#    result = func(*args)
#    return '%s says that %s%s = %s' % \
#        (current_process().name, func.__name__, args, result)
#
##
## Functions referenced by tasks
##
#
#def mul(a, b):
#    time.sleep(0.5*random.random())
#    return a * b
#
#def plus(a, b):
#    time.sleep(0.5*random.random())
#    return a + b
#
##
##
##
#
#def test():
#    NUMBER_OF_PROCESSES = 4
#    TASKS1 = [(mul, (i, 7)) for i in range(20)]
#    TASKS2 = [(plus, (i, 8)) for i in range(10)]
#
#    # Create queues
#    task_queue = Queue()
#    done_queue = Queue()
#
#    # Submit tasks
#    for task in TASKS1:
#        task_queue.put(task)
#
#    # Start worker processes
#    for i in range(NUMBER_OF_PROCESSES):
#        Process(target=worker, args=(task_queue, done_queue)).start()
#
#    # Get and print results
#    print 'Unordered results:'
#    for i in range(len(TASKS1)):
#        print '\t', done_queue.get()
#
#    # Add more tasks using `put()`
#    for task in TASKS2:
#        task_queue.put(task)
#
#    # Get and print some more results
#    for i in range(len(TASKS2)):
#        print '\t', done_queue.get()
#
#    # Tell child processes to stop
#    for i in range(NUMBER_OF_PROCESSES):
#        task_queue.put('STOP')
#
#
#if __name__ == '__main__':
#    freeze_support()
#    test()