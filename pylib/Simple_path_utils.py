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
## path_A       ==========     or   ============   or =====        (either are idx 0)
## path_B       ==========          ====              ===========



def path_A_contains_path_B(simple_path_A, simple_path_B):

    if len(simple_path_B) > len(simple_path_A):
        return False

    if simple_path_B[0] not in simple_path_A:
        return False
    
    idx_A = simple_path_A.index(simple_path_B[0])
    if idx_A < 0:
        return False

    if idx_A + len(simple_path_B) > len(simple_path_A):
        return False

    idx_B = 0
    while idx_B < len(simple_path_B):
        idx_B += 1
        idx_A += 1
        if (idx_B < len(simple_path_B) and
            idx_A < len(simple_path_A) and
            simple_path_A[idx_A] != simple_path_B[idx_B]):
            return False

    return True


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
                else: #scenario 5, return longer path
                    if len(simple_path_A) >= len(simple_path_B):
                        return simple_path_A
                    else:
                        return simple_path_B
            
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



def merge_adjacent_segments(segment_list):

    segment_list = sorted(segment_list, key=lambda x: x[0])
    
    ret_segments = list()
    ret_segments.append(list(segment_list.pop(0)))

    while len(segment_list) > 0:
        next_seg = list(segment_list.pop(0))

        prev_seg = ret_segments[-1]
        if next_seg[0] == prev_seg[1] + 1:
            # merge them
            assert(next_seg[1] > prev_seg[1])
            prev_seg[1] = next_seg[1]
        else:
            ret_segments.append(next_seg)

    return ret_segments




if __name__ == '__main__':

    ## run unit tests:

    # Tests that should return True

    path_a = ["n0", "n1", "n2", "n3", "n4", "n5", "n6"]
    path_b =             ["n2", "n3", "n4", "n5", "n6", "n7", "n8"]
    test = are_overlapping_and_compatible_NO_gaps_in_overlap(path_a, path_b)
    print("path_a: {}\npath_b: {}\ncompatible_NO_gaps: {}".format(path_a, path_b, test))
    assert(test is True)


    path_a = ["n0", "n1", "n2", "n3", "n4", "n5", "n6"]
    path_b =             ["n2", "n3", "n4", "n5"]
    test = are_overlapping_and_compatible_NO_gaps_in_overlap(path_a, path_b)
    print("path_a: {}\npath_b: {}\ncompatible_NO_gaps: {}".format(path_a, path_b, test))
    assert(test is True)

    
    path_a =              ["n2", "n3", "n4", "n5"]
    path_b = ["n0", "n1", "n2", "n3", "n4", "n5", "n6"]
    test = are_overlapping_and_compatible_NO_gaps_in_overlap(path_a, path_b)
    print("path_a: {}\npath_b: {}\ncompatible_NO_gaps: {}".format(path_a, path_b, test))
    assert(test is True)


    path_a =             ["n2", "n3", "n4", "n5", "n6", "n7", "n8"]
    path_b = ["n0", "n1", "n2", "n3", "n4", "n5", "n6"]
    test = are_overlapping_and_compatible_NO_gaps_in_overlap(path_a, path_b)
    print("path_a: {}\npath_b: {}\ncompatible_NO_gaps: {}".format(path_a, path_b, test))
    assert(test is True)

    

    ################################
    # Tests that should return False

    path_a =             ["n2", "n10", "n4", "n5", "n6", "n7", "n8"]
    path_b = ["n0", "n1", "n2", "n3", "n4", "n5", "n6"]
    test = are_overlapping_and_compatible_NO_gaps_in_overlap(path_a, path_b)
    print("path_a: {}\npath_b: {}\ncompatible_NO_gaps: {}".format(path_a, path_b, test))
    assert(test is False)
    
        
    path_a =             ["n2", "X10", "n4", "n5", "n6", "n7", "n8"]
    path_b = ["n0", "n1", "n2", "n3", "n4", "n5", "n6"]
    test = are_overlapping_and_compatible_NO_gaps_in_overlap(path_a, path_b)
    print("path_a: {}\npath_b: {}\ncompatible_NO_gaps: {}".format(path_a, path_b, test))
    assert(test is False)

    
    path_a = ["n2", "n10"]
    path_b =               ["n3", "n4", "n5", "n6"]
    test = are_overlapping_and_compatible_NO_gaps_in_overlap(path_a, path_b)
    print("path_a: {}\npath_b: {}\ncompatible_NO_gaps: {}".format(path_a, path_b, test))
    assert(test is False)
    

    ###################
    ## Test merging paths

    path_a = ["n1", "n2", "n3"]
    path_b =       ["n2", "n3", "n4", "n5"]
    merged_path = merge_simple_paths(path_a, path_b)
    print("path_a: {}\npath_b: {}\nmerged_path: {}".format(path_a, path_b, merged_path))
    assert(merged_path == ["n1", "n2", "n3", "n4", "n5"])



    path_a = ["n1", "n2", "n3"]
    path_b =       ["n2", "n3"]
    merged_path = merge_simple_paths(path_a, path_b)
    print("path_a: {}\npath_b: {}\nmerged_path: {}".format(path_a, path_b, merged_path))
    assert(merged_path == ["n1", "n2", "n3"])


    path_a =              ["n3", "n4", "n5"]
    path_b = ["n1", "n2", "n3"]
    merged_path = merge_simple_paths(path_a, path_b)
    print("path_a: {}\npath_b: {}\nmerged_path: {}".format(path_a, path_b, merged_path))
    assert(merged_path == ["n1", "n2", "n3", "n4", "n5"])


    
    ####################
    ## Test containments
    path_a = ["n1", "n2", "n3"]
    path_b =       ["n2"]
    test = path_A_contains_path_B(path_a, path_b)
    print("path_a: {} contains path_b: {} = {}".format(path_a, path_b, test))
    assert(test is True)


    path_a = ["n1", "n2", "n3"]
    path_b = ["n1", "n2", "n3"]
    test = path_A_contains_path_B(path_a, path_b)
    print("path_a: {} contains path_b: {} = {}".format(path_a, path_b, test))
    assert(test is True)


    path_a = ["n1", "n2", "n3"]
    path_b = ["n1", "n2", "n3", "n4"]
    test = path_A_contains_path_B(path_a, path_b)
    print("path_a: {} contains path_b: {} = {}".format(path_a, path_b, test))
    assert(test is False)

    path_a = ["n1", "n2", "n3"]
    path_b = ["n0", "n1", "n2", "n3", "n4"]
    test = path_A_contains_path_B(path_a, path_b)
    print("path_a: {} contains path_b: {} = {}".format(path_a, path_b, test))
    assert(test is False)


    path_a = ["n1", "n2", "n3"]
    path_b = ["n3", "n4"]
    test = path_A_contains_path_B(path_a, path_b)
    print("path_a: {} contains path_b: {} = {}".format(path_a, path_b, test))
    assert(test is False)


        
    sys.exit(0)

    
