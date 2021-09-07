#!/usr/bin/env python3

import time
import multiprocessing
import random
import logging


logger = logging.getLogger(__name__)

SLEEPTIME = 0.1

class MultiProcessManager:

    def __init__(self, num_parallel_processes, queue=None):
        
        self.num_parallel_processes = num_parallel_processes
        self.num_running = 0
        self.num_successes = 0
        self.num_errors = 0
        self.process_list = list()
        self.queue = queue
        self.captured_queue_contents = list()
        
    def launch_process(self, process):

        logger.debug("-launching process")
        
        if self.num_running >= self.num_parallel_processes:
            self.wait_for_open_slot()

        process.start()
        self.process_list.append(process)
        self.num_running += 1


    def wait_for_open_slot(self):
        logger.debug("-waiting for open slot")

        counter = 0
        while self.num_running >= self.num_parallel_processes:
            logger.debug("\twaiting for open slot round({})".format(counter))
            counter += 1
            self._screen_running_processes()
            time.sleep(SLEEPTIME)

    def _screen_running_processes(self):

        logger.debug("-screening running processes")

        if self.queue is not None and not self.queue.empty():
            logger.debug("\t-reaping queue")
            while not self.queue.empty():
                logger.debug("\t\t-try reaping entry from queue")
                entry = self.queue.get()
                logger.debug("\t\t-reaped entry")
                self.captured_queue_contents.append(entry)
                
            logger.debug("\t\t-done reaping queue")

        
        completed_processes = list()
            
        for i, process in enumerate(self.process_list):
            if process.is_alive():
                logger.debug("\t-process {} is alive.".format(i))
            else:
                logger.debug("\t-process {} is finished.".format(i))
                completed_processes.append(i)

        if completed_processes:
            logger.debug("\t-processing {} completed processes.".format(len(completed_processes)))
            completed_processes.reverse()
            for completed_process_idx in completed_processes:
                
                process = self.process_list[completed_process_idx]
                logger.debug("\t\t\t<joining process {}>".format(completed_process_idx))
                process.join()
                logger.debug("\t\t\t<joined process {} having exit code {}>".format(completed_process_idx, process.exitcode))
                if process.exitcode == 0:
                    self.num_successes += 1
                else:
                    self.num_errors += 1
                    logger.debug("-captured a failed process")

                del self.process_list[completed_process_idx]
                self.num_running -= 1

                
        logger.debug("\t-done screening running processes.")
        
                
    def wait_for_remaining_processes(self):

        logger.debug("-waiting for remaining processes")
        
        while self.num_running > 0:
            logger.debug("-waiting on {} processes".format(self.num_running))
            self._screen_running_processes()
            time.sleep(SLEEPTIME)

        logger.debug("-done waiting. All processes are completed")

        return self.num_errors
        

    def summarize_status(self):
        return("{} jobs succeeded & {} jobs failed".format(self.num_successes, self.num_errors))


    def retrieve_queue_contents(self):
        return self.captured_queue_contents


def set_debug():
    logger.setLevel(logging.DEBUG)

        

def test_mpm(num_parallel_processes=8, num_total_processes=100):

    def runner(id,q):
        print("running id:{}".format(id))
        x = id / (id % 10)  # should error as div-by-zero on occasion
        time.sleep(random.randint(0,10))
        q.put(id)
        
    q = multiprocessing.Queue()
    mpm = MultiProcessManager(num_parallel_processes, q)
        
    for i in range(num_total_processes):
    
        p = multiprocessing.Process(target=runner, args=(i,q))

        mpm.launch_process(p)


    mpm.wait_for_remaining_processes()
    
    



if __name__=='__main__':

    # run test
    logging.basicConfig(level=logging.INFO, 
                        format='%(asctime)s : %(levelname)s : %(message)s',
                        datefmt='%H:%M:%S')

    logger.setLevel(logging.DEBUG)
    test_mpm()
    
