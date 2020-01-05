import sys, os, re
import Splice_graph
import MultiPath
import MultiPathCounter
import networkx as nx
import Simple_path_utils as Simple_path_utils

class MultiPathGraphNode:

    mp_id_counter = 0
    
    def __init__(self, multiPath, count, lend, rend, mpg):
        self._multiPath = multiPath
        self._count = count
        self._lend = lend
        self._rend = rend

        self._mpg = mpg

        self._containments = set() # other MPG nodes fully contained by this node.

        MultiPathGraphNode.mp_id_counter += 1
        self._id = "mp{}".format(MultiPathGraphNode.mp_id_counter)


    def get_multiPath(self):
        return self._multiPath

    def get_multipath_graph(self):
        return self._mpg

    def get_simple_path(self):
        return(self._multiPath.get_path())

    def get_coords(self):
        return(self._lend, self._rend)

    def get_count(self):
        return self._count

    
    def __repr__(self):
        return("mp:{} {}-{} C:{}".format(self.get_simple_path(), self._lend, self._rend, self._count))
    

    def coords_overlap(self, other_node):
        my_lend, my_rend = self.get_coords()
        other_lend, other_rend = other_node.get_coords()

        if my_lend <= other_rend and my_rend >= other_lend:
            return True
        else:
            return False

        
    def add_containment(self, other_node):
        self._containments.add(other_node)

    
    def get_containments(self):
        return(list(self._containments))

    
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

        
    def get_mpgn_pair_token(mpgn_A, mpgn_B):
        assert(type(mpgn_A) == MultiPathGraphNode)
        assert(type(mpgn_B) == MultiPathGraphNode)
        
        token = "^^".join(sorted((mpgn_A._id, mpgn_B._id)))
        return token
        
        
class MultiPathGraph:


    def __init__(self, multiPathCounter, splice_graph):

        assert(type(multiPathCounter) == MultiPathCounter.MultiPathCounter)
        assert(type(splice_graph) == Splice_graph.Splice_graph)
        
        self._splice_graph = splice_graph
        
        mp_graph = nx.DiGraph()
        self._mp_graph = mp_graph

        self._mp_graph_nodes_list = list()

        self._incompatible_pairs = set() # store sets of pairs of incompatible nodes.
        
        multiPathCountPairs = multiPathCounter.get_all_MultiPathCountPairs()
        for mpCountPair in multiPathCountPairs:
            mp, count = mpCountPair.get_multipath_and_count()

            first_node_id = mp[0]
            last_node_id = mp[-1]

            first_node_obj = splice_graph.get_node_obj_via_id(first_node_id)
            last_node_obj = splice_graph.get_node_obj_via_id(last_node_id)

            lend_coord = first_node_obj.get_coords()[0]
            rend_coord = last_node_obj.get_coords()[1]

            mp_graph_node = MultiPathGraphNode(mp, count, lend_coord, rend_coord, mpg=self)
            mp_graph.add_node(mp_graph_node)

            self._mp_graph_nodes_list.append(mp_graph_node)

            
        ## sort by rend
        self._mp_graph_nodes_list = sorted(self._mp_graph_nodes_list, key=lambda x: x._rend)

        ordered_nodes = self._mp_graph_nodes_list
        
        ## define edges, containments, and incompatibilities
        for i in range(0, len(ordered_nodes)):
            node_i = ordered_nodes[i]

            for j in range(i-1, -1, -1):
                node_j = ordered_nodes[j]

                #print("comparing {},{}".format(node_j, node_i))
                
                if node_j._rend < node_i._lend:
                    break # all earlier node j's will also be non-overlapping

                if node_i.contains_other_node(node_j):
                    node_i.add_containment(node_j)
                elif node_j.contains_other_node(node_i):
                    node_j.add_containment(node_i)
                elif node_i.compatible(node_j):
                    self._mp_graph.add_edge(node_j, node_i)
                    print("adding edge: {},{}".format(node_j, node_i))
                else:
                    # incompatible pairs
                    incompatible_pair_token = MultiPathGraphNode.get_mpgn_pair_token(node_i, node_j)
                    self._incompatible_pairs.add(incompatible_pair_token)

    
    def get_ordered_nodes(self):
        return list(self._mp_graph_nodes_list)


    
    def has_edge(self, multiPath_before, multiPath_after):
        return self._mp_graph.has_edge(multiPath_before, multiPath_after)

    
    def incompatible_mpgn_pair(self, mpgn_A, mpgn_B):
        mpgn_pair_token = MultiPathGraphNode.get_mpgn_pair_token(mpgn_A, mpgn_B)
        if mpgn_pair_token in self._incompatible_pairs:
            return True
        else:
            return False

        
