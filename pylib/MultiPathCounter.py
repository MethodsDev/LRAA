#!/usr/bin/env python
# encoding: utf-8

import sys, os, re
from collections import defaultdict


class MultiPathCountPair:

    def __init__(self, multipath, count=1):
        self._multipath = multipath
        self._count = count
        return

    def increment(self, increment=1):
        self._count += increment
        return

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
            self._multipath_counter[multipath_key].increment()
        else:
            self._multipath_counter[multipath_key] = MultiPathCountPair(multipath)
        
        return

    
    def __repr__(self):

        ret_text = "\t"

        for multipath_count_pair in self._multipath_counter.values():
            ret_text += str(multipath_count_pair) + "\n"

        return ret_text
            

    
