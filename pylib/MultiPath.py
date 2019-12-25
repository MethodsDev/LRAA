#!/usr/bin/env python
# encoding: utf-8

import sys, os, re
from collections import defaultdict
from PASA_SALRAA_Globals import SPACER
import Simple_path_utils

class MultiPath:
    
    
    def __init__(self, splice_graph, paths_list):

        self._splice_graph = splice_graph

        self._multipath = self._merge_paths_to_multi_path(paths_list)

        return


    def __len__(self):
        return(len(self._multipath))


    def __getitem__(self, i):
        length = len(self)
        if i < 0:
            i += length
        if 0 <= i < length:
            return self._multipath[i]
        
        raise IndexError('Index out of range: {}'.format(i))
    
    def _merge_paths_to_multi_path(self, paths_list):


        paths_to_asm = [path for path in paths_list] # copy incoming list

        ## perform cycles of merge attempts, retaining relative ordering.
        assembled_paths = list()        
        # seed it with the first entry
        assembled_paths.append(paths_to_asm.pop(0))
    
        
        while True:
            
            unmerged_paths = list()

            seed_asm = assembled_paths[-1]

            merged_flag = False
            
            for other_path in paths_to_asm:
                # check for merge
                if Simple_path_utils.are_overlapping_and_compatible_NO_gaps_in_overlap(seed_asm, other_path):
                    #mergeable, so lets merge them
                    merged_asm = Simple_path_utils.merge_simple_paths(seed_asm, other_path)
                    # update seed asm
                    seed_asm = assembled_paths[-1] = merged_asm
                    merged_flag = True
                else:
                    unmerged_paths.append(other_path)
                    
            if unmerged_paths:
                if not merged_flag:
                    # must make a new seed.
                    assembled_paths.append(unmerged_paths.pop(0))

                paths_to_asm = unmerged_paths

            else:
                break # done assembling
        
                    
        ## build multipath for post-assembly products
        multipath = []
        for i, path in enumerate(assembled_paths):
            multipath += path
            if i != len(assembled_paths) -1:
                multipath.append(SPACER)

        return multipath


    def __repr__(self):
        return(str(self._multipath))
    

    

if __name__ == '__main__':

    
    
    ## unit tests

    paths_list = [  ["n1", "n2"],
                          ["n2", "n3"]  ]
    mp = MultiPath._merge_paths_to_multi_path(None, paths_list)
    print("Path list: {} merged into multipath: {}".format(paths_list, mp))
    assert(mp == ["n1", "n2", "n3"])


    paths_list = [  ["n1", "n2"],
                          ["n3", "n4"]  ]
    mp = MultiPath._merge_paths_to_multi_path(None, paths_list)
    print("Path list: {} merged into multipath: {}".format(paths_list, mp))
    assert(mp == ['n1', 'n2', '???', 'n3', 'n4'] )

       
    sys.exit(0)

    
