#!/usr/bin/env python
# encoding: utf-8

import sys, os, re
from collections import defaultdict
import logging

logger = logging.getLogger(__name__)

class MultiPathCountPair:

    def __init__(self, multipath, count=1):
        self._multipath = multipath
        self._count = count
        return

    def increment(self, increment=1):
        self._count += increment
        return

    def get_multipath_and_count(self):
        return(self._multipath, self._count)

    def include_read_type(self, read_type):
        self._multipath.include_read_type(read_type)

    def include_read_name(self, read_name):
        self._multipath.include_read_name(read_name)
    
        
    def __repr__(self):
        ret_text = "{}\t{}".format(str(self._multipath), self._count)
        return ret_text
    
    
class MultiPathCounter:
    
    
    def __init__(self):

        self._multipath_counter = dict()
        
        return


    def add(self, multipath):

        multipath_key = str(multipath)
        if multipath_key in self._multipath_counter:
            mp = self._multipath_counter[multipath_key]
        else:
            mp = self._multipath_counter[multipath_key] = MultiPathCountPair(multipath)

        
        mp.increment()
        mp.include_read_type(multipath.get_read_types())
        mp.include_read_name(multipath.get_read_names())
                    
        return
    

    
    def get_all_MultiPathCountPairs(self):
        return(self._multipath_counter.values())

    
    def __repr__(self):

        ret_text = "\t"

        for multipath_count_pair in self._multipath_counter.values():
            ret_text += str(multipath_count_pair) + "\n"

        return ret_text
            

    
