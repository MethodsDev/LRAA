#!/usr/bin/env python
# encoding: utf-8

import sys, os, re
from collections import defaultdict
from PASA_SALRAA_Globals import SPACER

# Namespace: Simple_path_utils
# includes basic functions for evaluating relationships between simple paths in the graph.



## Scenarios:
#
## Scenario 1
## path_A  ==============
## path_B        ===============  (idx_B == 0)
## or
#
## Scenario 2
## path_A        ================= (idx_A == 0)
## path_B  ==============
## or
#
## Scenarion 3
## path_A  ======================= (idx_B == 0)
## path_B         =======
## or
#
#  Scenario 4
## path_A         =======
## path_B  ======================= (idx_A == 0)
## or
#
# Scenario 5
## path_A       ==========     (either are idx 0)
## path_B       ==========



def are_overlapping_and_compatible_NO_gaps_in_overlap(simple_path_A, simple_path_B):
    
    ## find first non-spacer match between two paths.  Ensure remaining parts of paths are identical
    for idx_A in range(0, len(simple_path_A)):
        A_node_id = simple_path_A[idx_A]
        if A_node_id == SPACER:
            continue
        if A_node_id in simple_path_B:
            idx_B = simple_path_B.index(A_node_id)

            ## one of the indexes needs to start at zero or there'll be some unmatched upstream nodes.
            if (idx_B != 0 and idx_A != 0):
                return False

            # ensure remainder of paths overlap, no gaps allowed.
            idx_A += 1
            idx_B += 1
            while (idx_A < len(simple_path_A) and idx_B < len(simple_path_B)):
                if simple_path_A[idx_A] == SPACER or simple_path_A[idx_A] != simple_path_B[idx_B]:
                    return False
                idx_A += 1
                idx_B += 1
                
            # if made it here, the remainder of the paths are gap-free and identical
            return True

    return False # no matching index between paths




def are_overlapping_and_compatibile_ALLOW_gaps_in_overlap(simple_path_A, simple_path_B):
    pass




def merge_simple_paths(simple_path_A, simple_path_B):
    
    if not are_overlapping_and_compatible_NO_gaps_in_overlap(simple_path_A, simple_path_B):
        raise RuntimeException("cannot merge paths that are not compatible in overlapping region")


    ## find first non-spacer match between two paths, then merge.
    for idx_A in range(0, len(simple_path_A)):
        A_node_id = simple_path_A[idx_A]
        if A_node_id == SPACER:
            continue
        if A_node_id in simple_path_B:
            idx_B = simple_path_B.index(A_node_id)

            merged_path = None
            if idx_A == 0: # scenarios 2,4,5
                if idx_B > 0: # scenario 2 or 4
                    # begin path with B prefix
                    merged_path = simple_path_B
                    # if A extends past B, need to include that.
                    #  path A:        0 1 2 3 4 5 6
                    #  path B:    0 1 2 3 4 5 6
                    extension_idx_A = len(simple_path_B) - idx_B
                    if extension_idx_A < len(simple_path_A):
                        merged_path.extend(simple_path_A[extension_idx_A: ])
                    return merged_path
                else: #scenario 5, must be identical
                    if simple_path_A == simple_path_B:
                        return simple_path_A
                    else:
                        raise RuntimeError("Error, paths are not identical: {} and {}".format(simple_path_A, simple_path_B))
            
            else:
                # idx_A != 0, so idx_B must be == 0
                assert(idx_B == 0)
                # scenarios 1,3
                # begin path with A prefix
                merged_path = simple_path_A
                # if A extends past B, need to include that.
                #  path A:   0 1 2 3 4 5 6
                #  path B:         0 1 2 3 4 5 6
                extension_idx_B = len(simple_path_A) - idx_A
                if extension_idx_B < len(simple_path_B):
                    merged_path.extend(simple_path_B[extension_idx_B: ])
                return merged_path
            
    raise RuntimeException("Error, could not merge simple paths {} and {} ... bug... ".format(simple_path_A, simple_path_B))





if __name__ == '__main__':

    ## run unit tests:
    path_a = ["n0", "n1", "n2", "n3", "n4", "n5", "n6"]
    path_b =             ["n2", "n3", "n4", "n5", "n6", "n7", "n8"]

    test = are_overlapping_and_compatible_NO_gaps_in_overlap(path_a, path_b)
    print("path_a: {}\npath_b: {}\ncompatible_NO_gaps: {}".format(path_a, path_b, test))
    
    sys.exit(0)

    
