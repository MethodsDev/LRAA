import sys, os, re
import Splice_graph
import MultiPath
import MultiPathGraph
import MultiPathCounter
import networkx as nx
import Simple_path_utils as Simple_path_utils
from PASA_SALRAA_Globals import SPACER
from GenomeFeature import Intron, Exon

import logging

logger = logging.getLogger(__name__)


class MultiPathGraphNode:

    mp_id_counter = 0
    
    def __init__(self, multiPathObj, count, lend, rend, mpg):
        self._multiPath = multiPathObj
        self._count = count
        self._weight = 1.0 # weight applied to read count, varied in pasa-salraa
        
        self._lend = lend
        self._rend = rend

        self._mpg = mpg

        self._containments = set() # other MPG nodes fully contained by this node.

        MultiPathGraphNode.mp_id_counter += 1
        self._id = "mp{}".format(MultiPathGraphNode.mp_id_counter)

        self._seq_length = self._compute_seq_length()
        


    def get_id(self):
        return self._id
    
    def get_multiPathObj(self):
        return self._multiPath

    def get_multipath_graph(self):
        return self._mpg

    def get_simple_path(self):
        return(self._multiPath.get_simple_path())

    def get_coords(self):
        return(self._lend, self._rend)

    def get_count(self):
        return self._count

    def set_count(self, count_val):
        self._count = count_val
        return

    def get_weight(self):
        return self._weight

    def set_weight(self, weight):
        self._weight = weight
        return

    def get_score(self):
        # score reflects density of read evidence in this mpgn 
        return ( (self._count * self._weight) / self._seq_length )
    
    def get_seq_length(self):
        return self._seq_length

    
    def __repr__(self):
        return("mp:{} {}-{} C:{} len:{}".format(self.get_simple_path(), self._lend, self._rend, self._count, self._seq_length))

        
    def has_successors(self):
        if len(list(self._mpg.successors(self))) > 0:
            return True
        else:
            return False

    def get_successors(self):
        return(list(self._mpg.successors(self)))

    def has_predecessors(self):
        if len(list(self._mpg.predecessors(self))) > 0:
            return True
        else:
            return False

    def get_predecessors(self):
        return(list(self._mpg.predecessors(self)))

    
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


    def get_splice_graph(self):
        return self.get_multipath_graph().get_splice_graph()


    def get_splice_graph_node_objs_for_path(self):
        simple_path = self.get_simple_path()
        splice_graph_node_objs = list()

        sg = self.get_splice_graph()
        
        for node_id in simple_path:
            if node_id == SPACER:
                splice_graph_node_objs.append(None)
            else:
                splice_graph_node_obj = sg.get_node_obj_via_id(node_id)
                splice_graph_node_objs.append(splice_graph_node_obj)

        return splice_graph_node_objs


    def _compute_seq_length(self):
        sg_nodes = self.get_splice_graph_node_objs_for_path()

        seq_length = 0
        for sg_node in sg_nodes:
            if sg_node is not None and type(sg_node) == Exon:
                lend, rend = sg_node.get_coords()
                exon_len = rend - lend + 1
                seq_length += exon_len

        return seq_length
        
