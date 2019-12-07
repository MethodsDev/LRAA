#!/usr/bin/env python
# encoding: utf-8

import sys, os, re
from collections import defaultdict


class MultiPath:
    
    
    def __init__(self, splice_graph, paths_list):

        self._splice_graph = splice_graph

        self._multipath = self._merge_paths_to_multi_path(paths_list)

        return

    def _merge_paths_to_multi_path(self, paths_list):

        ## //FIXME: proper implementation needed.  Using placeholder for now.
        multipath = []
        for path in paths_list:
            multipath += path + ["???"]

        return multipath


    def __repr__(self):
        return(str(self._multipath))
    

    
