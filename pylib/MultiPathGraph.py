import sys, os, re
import Splice_graph
import MultiPath
import MultiPathCounter
import networkx as nx
import Simple_path_utils as Simple_path_utils
from PASA_SALRAA_Globals import SPACER
from GenomeFeature import Intron, Exon
from MultiPathGraphNode import MultiPathGraphNode

import logging

logger = logging.getLogger(__name__)

    
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
                    logger.debug("adding edge: {},{}".format(node_j, node_i))
                else:
                    # incompatible pairs
                    incompatible_pair_token = MultiPathGraphNode.get_mpgn_pair_token(node_i, node_j)
                    self._incompatible_pairs.add(incompatible_pair_token)

    
    def get_ordered_nodes(self):
        # these are sorted by rend
        return list(self._mp_graph_nodes_list)

    
    def has_edge(self, multiPath_before, multiPath_after):
        return self._mp_graph.has_edge(multiPath_before, multiPath_after)


    def successors(self, mpgn):
        return self._mp_graph.successors(mpgn)

    def predecessors(self, mpgn):
        return self._mp_graph.predecessors(mpgn)

    
    def incompatible_mpgn_pair(self, mpgn_A, mpgn_B):
        mpgn_pair_token = MultiPathGraphNode.get_mpgn_pair_token(mpgn_A, mpgn_B)
        if mpgn_pair_token in self._incompatible_pairs:
            return True
        else:
            return False

        
    def get_splice_graph(self):
        return self._splice_graph


    def define_disjoint_graph_components(self):

        mpgn_list = self.get_ordered_nodes()
        mpgn_seen = set()

        component_list = list()
        
        while len(mpgn_list) > 0:

            queue = list()
            seed_node = mpgn_list.pop(0)

            if seed_node in mpgn_seen:
                continue

            # start a new component.
            component = list()
            queue.append(seed_node)

            while len(queue) > 0:
                node = queue.pop(0)
                if node not in mpgn_seen:
                    mpgn_seen.add(node)
                    component.append(node)
                    # add predecessors and successors to queue
                    if node.has_predecessors():
                        queue.extend(node.get_predecessors())
                    if node.has_successors():
                        queue.extend(node.get_successors())

            if len(component) > 0:
                component_list.append(component)

        logger.info("identified {} disjoint graph components".format(len(component_list)))

        return component_list

    
