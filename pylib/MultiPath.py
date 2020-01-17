#!/usr/bin/env python
# encoding: utf-8

import sys, os, re
from collections import defaultdict
from PASA_SALRAA_Globals import SPACER
import Simple_path_utils
from Util_funcs import coordpairs_overlap
import logging
from GenomeFeature import Exon

logger = logging.getLogger(__name__)


##
## The MultiPath is a list of node IDs that can be interrupted by SPACERs.
## The MultiPath stores the splice graph so that the nodes corresponding to those IDs can be retrieved as needed.
##


class MultiPath:
    
    
    def __init__(self, splice_graph, paths_list):

        self._splice_graph = splice_graph

        self._simple_path = self._merge_paths_to_simple_multi_path(paths_list)

        # determine span
        coords = list()
        for node_id in self._simple_path:
            if node_id != SPACER:
                coords.extend( self._splice_graph.get_node_obj_via_id(node_id).get_coords())
        coords = sorted(coords)

        self._lend = coords[0]
        self._rend = coords[-1]
        
        return


    def get_simple_path(self):
        return(list(self._simple_path)) # send a copy

    def get_splice_graph(self):
        return self._splice_graph
    
    def get_ordered_exons_and_introns(self):
        simple_path = self.get_simple_path()

        sg = self.get_splice_graph()
        
        # spacers become None
        
        exons_and_introns = list()
        for node_id in simple_path:
            if node_id == SPACER:
                exons_and_introns.append(None)
            else:
                obj = sg.get_node_obj_via_id(node_id)
                exons_and_introns.append(obj)

        return exons_and_introns

    
    def __len__(self):
        return(len(self._simple_path))


    def __getitem__(self, i):
        length = len(self)
        if i < 0:
            i += length
        if 0 <= i < length:
            return self._simple_path[i]
        
        raise IndexError('Index out of range: {}'.format(i))

    
    def get_coords(self):
        return(self._lend, self._rend)


    def get_cdna_length(self):

        cdna_len = 0
        
        exons_and_introns = self.get_ordered_exons_and_introns()
        for feature in exons_and_introns:
            if type(feature) == Exon:
                lend,rend = feature.get_coords()
                exon_len = rend - lend + 1
                cdna_len += exon_len

        return cdna_len
    
    
    
    def _merge_paths_to_simple_multi_path(self, paths_list):


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
        simple_multipath = []
        for i, path in enumerate(assembled_paths):
            simple_multipath += path
            if i != len(assembled_paths) -1:
                simple_multipath.append(SPACER)

        return simple_multipath



    def is_overlapping_contained_and_compatible(self, other_multipath):

        my_lend, my_rend = self.get_coords()
        other_lend, other_rend = other_multipath.get_coords()
        if (my_lend <= other_lend and my_rend > other_rend):
            return self.is_overlapping_and_compatible(other_multipath)
        else:
            return False

        

    def is_overlapping_and_compatible(self, other_multipath):
        # overlapping parts are required to be identical
        # compatible means no conflict detected.

        assert(type(other_multipath) == MultiPath)
        
        if not coordpairs_overlap(self.get_coords(), other_multipath.get_coords()):
            return False

        my_path = self.get_simple_path()
        other_path = other_multipath.get_simple_path()

        while(len(my_path) > 0 and len(other_path) > 0):
            my_node_id = my_path[0]
            if my_node_id == SPACER:
                my_path.pop(0)
                continue
            
            other_node_id = other_path[0]
            if other_node_id == SPACER:
                other_path.pop(0)
                continue

            if my_node_id == other_node_id:
                # great! advance and continue
                my_path.pop(0)
                other_path.pop(0)
                continue
            
            my_node_obj = self._splice_graph.get_node_obj_via_id(my_node_id)
            other_node_obj = self._splice_graph.get_node_obj_via_id(other_node_id)
            
            if coordpairs_overlap(my_node_obj.get_coords(), other_node_obj.get_coords()):
                # uh oh, they overlap but they're not the same.
                return False

            # advance the side that's before the other
            if my_node_obj.get_coords()[0] < other_node_obj.get_coords()[0]:
                my_path.pop(0)
            else:
                other_path.pop(0)

        return True # no conflicts detected
    

    def __repr__(self):
        return(str(self._simple_path))
    

    ################
    ## Class methods
    ################
        
    def multiPath_from_mpgn_list(mpgn_list):

        assert(len(mpgn_list) > 0)
        
        path_list = list()

        sg = mpgn_list[0].get_splice_graph()

        for mpgn in mpgn_list:
            path_list.append(mpgn.get_simple_path())

        mp_obj = MultiPath(sg, path_list)

        return mp_obj

    
        


        

if __name__ == '__main__':

    
    
    ## unit tests

    paths_list = [  ["n1", "n2"],
                          ["n2", "n3"]  ]
    mp = MultiPath._merge_paths_to_simple_multi_path(None, paths_list)
    print("Path list: {} merged into multipath: {}".format(paths_list, mp))
    assert(mp == ["n1", "n2", "n3"])


    paths_list = [  ["n1", "n2"],
                          ["n3", "n4"]  ]
    mp = MultiPath._merge_paths_to_simple_multi_path(None, paths_list)
    print("Path list: {} merged into multipath: {}".format(paths_list, mp))
    assert(mp == ['n1', 'n2', '???', 'n3', 'n4'] )

       
    sys.exit(0)

    
