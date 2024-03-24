#!/usr/bin/env python3

import time
import multiprocessing
import random
import logging

logger = logging.getLogger(__name__)

SLEEPTIME = 0.01

MPM_DEBUG = False

class MultiProcessManager:
    
    def __init__(self, num_parallel_processes, queue=None):
        
        self.num_parallel_processes = num_parallel_processes
        self.num_running = 0
        self.num_successes = 0
        self.num_errors = 0
        self.process_list = set()
        self.queue = queue
        self.captured_queue_contents = list()
        self.process_name_to_start_time = dict()
        
    def launch_process(self, process):

        if MPM_DEBUG:
            logger.debug("-launching process")
        
        if self.num_running >= self.num_parallel_processes:
            self.wait_for_open_slot()

        process.start()
        self.process_name_to_start_time[process.name] = time.time()
        self.process_list.add(process)
        self.num_running += 1


    def wait_for_open_slot(self):
        if MPM_DEBUG:
            logger.debug("-waiting for open slot")

        counter = 0
        while self.num_running >= self.num_parallel_processes:
            if MPM_DEBUG:
                logger.debug("\twaiting for open slot round({})".format(counter))
            counter += 1
            self._screen_running_processes()
            time.sleep(SLEEPTIME)

    def _screen_running_processes(self):

        logger.debug("-screening running processes")

        if self.queue is not None and not self.queue.empty():
            if MPM_DEBUG:
                logger.debug("\t-reaping queue")
            while not self.queue.empty():
                if MPM_DEBUG:
                    logger.debug("\t\t-try reaping entry from queue")

                entry = self.queue.get()

                if MPM_DEBUG:
                    logger.debug("\t\t-reaped entry")

                self.captured_queue_contents.append(entry)

            if MPM_DEBUG:
                logger.debug("\t\t-done reaping queue")

        
        completed_processes = set()
            
        for process in self.process_list:
            if process.is_alive():
                if MPM_DEBUG:
                    logger.debug("\t-process {} is alive.".format(process.name))
            else:
                if MPM_DEBUG:
                    logger.debug("\t-process {} is finished.".format(process.name))
                completed_processes.add(process)

        if completed_processes:
            if MPM_DEBUG:
                logger.debug("\t-processing {} completed processes.".format(len(completed_processes)))

            for completed_process in completed_processes:
                
                if MPM_DEBUG:
                    logger.debug("\t\t\t<joining process {}>".format(completed_process.name))
                completed_process.join()
                if MPM_DEBUG:
                    logger.debug("\t\t\t<joined process {} having exit code {}>".format(completed_process.name, completed_process.exitcode))
                if completed_process.exitcode == 0:
                    self.num_successes += 1
                else:
                    self.num_errors += 1
                    if MPM_DEBUG:
                        logger.debug("-captured a failed process")

                # remove completed process from process list
                self.process_list.remove(completed_process)
                self.num_running -= 1

                
        logger.debug("\t-done screening running processes.")
        
                
    def wait_for_remaining_processes(self):

        if MPM_DEBUG:
            logger.debug("-waiting for remaining processes")
        
        while self.num_running > 0:
            if MPM_DEBUG:
                logger.debug("-waiting on {} processes".format(self.num_running))
            self._screen_running_processes()
            time.sleep(SLEEPTIME)

        if MPM_DEBUG:
            logger.debug("-done waiting. All processes are completed")

        return self.num_errors
        

    def summarize_status(self):
        return("{} jobs succeeded & {} jobs failed".format(self.num_successes, self.num_errors))


    def retrieve_queue_contents(self):
        return self.captured_queue_contents


def set_debug():
    logger.setLevel(logging.DEBUG)
    global MPM_DEBUG
    MPM_DEBUG = True


def test_mpm(num_parallel_processes=8, num_total_processes=100):

    def runner(id,q):
        logger.info("running id:{}".format(id))
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
    
