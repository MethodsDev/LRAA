import sys, os, re
import Splice_graph
import MultiPath
import MultiPathCounter
import networkx as nx
import Simple_path_utils as Simple_path_utils

class MultiPathGraphNode:

    mp_id_counter = 0
    
    def __init__(self, multiPath, count, lend, rend):
        self._multiPath = multiPath
        self._count = count
        self._lend = lend
        self._rend = rend

        self._containments = set() # other MPG nodes fully contained by this node.

        MultiPathGraphNode.mp_id_counter += 1
        self._id = "mp{}".format(MultiPathGraphNode.mp_id_counter)
        

    def get_simple_path(self):
        return(self._multiPath.get_path())

    def get_coords(self):
        return(self._lend, self._rend)


    def coords_overlap(self, other_node):
        my_lend, my_rend = self.get_coords()
        other_lend, other_rend = other_node.get_coords()

        if my_lend <= other_rend and my_rend >= other_lend:
            return True
        else:
            return False

        
    def add_containment(self, other_node):
        self._containments.add(other_node)

        
    def contains_other_node(self, other_node):
        if Simple_path_utils.path_A_contains_path_B(self.get_simple_path(), other_node.get_simple_path()):
            return True
        else:
            return False

    def compatible(self, other_node):
        if self._multiPath.is_overlapping_and_compatible(other_node._multiPath):
            return True
        else:
            return False

        
class MultiPathGraph:


    def __init__(self, multiPathCounter, splice_graph):

        self._splice_graph = splice_graph
        
        mp_graph = nx.DiGraph()
        self._mp_graph = mp_graph

        self._mp_graph_nodes_list = list()
        
        multiPathCountPairs = multiPathCounter.get_all_MultiPathCountPairs()
        for mpCountPair in multiPathCountPairs:
            mp, count = mpCountPair.get_multipath_and_count()

            first_node_id = mp[0]
            last_node_id = mp[-1]

            first_node_obj = splice_graph.get_node_obj_via_id(first_node_id)
            last_node_obj = splice_graph.get_node_obj_via_id(last_node_id)

            lend_coord = first_node_obj.get_coords()[0]
            rend_coord = last_node_obj.get_coords()[1]

            mp_graph_node = MultiPathGraphNode(mp, count, lend_coord, rend_coord)
            mp_graph.add_node(mp_graph_node)

            self._mp_graph_nodes_list.append(mp_graph_node)


        self._mp_graph_nodes_list = sorted(self._mp_graph_nodes_list, key=lambda x: x._rend)

        ordered_nodes = self._mp_graph_nodes_list
        
        ## define edges, containments, and incompatibilities
        for i in range(0, len(ordered_nodes)):
            node_i = ordered_nodes[i]
            for j in range(i, -1, -1):
                node_j = ordered_nodes[j]

                if node_j._rend < node_i._lend:
                    break # all earlier node j's will also be non-overlapping

                if node_i.contains_other_node(node_j):
                    node_i.add_containment(node_j)
                elif node_i.compatible(node_j):
                    self._mp_graph.add_edge(node_j, node_i)
                    
            
        
