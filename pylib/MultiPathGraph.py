import sys, os, re
import Splice_graph
import MultiPath
import MultiPathCounter
import networkx as nx


class MultiPathGraphNode:

    def __init__(self, multiPath, count, lend, rend):
        self._multiPath = multiPath
        self._count = count
        self._lend = lend
        self._rend = rend

        self._containments = list() # other MPG nodes fully contained by this node.


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


        
            
