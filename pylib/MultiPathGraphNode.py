import sys, os, re
import Splice_graph
import MultiPath
import MultiPathGraph
import MultiPathCounter
import networkx as nx
import Simple_path_utils as Simple_path_utils
from PASA_SALRAA_Globals import SPACER
from GenomeFeature import Intron, Exon
import math

import logging

logger = logging.getLogger(__name__)


NORMALIZE_SCORES_BY_SEQ_LENGTH = False
MIN_WEIGHT = 0.01

class MultiPathGraphNode:

    mp_id_counter = 0
    
    def __init__(self, multiPathObj, count, lend, rend, mpg):
        self._multiPath = multiPathObj
        self._count = count
        self._weight = 1.0 # weight applied to read count, varied in pasa-salraa
        self._prev_weight = 1.0 # retain previous weight setting
        
        self._lend = lend
        self._rend = rend

        self._mpg = mpg

        self._containments = set() # other MPG nodes fully contained by this node.

        MultiPathGraphNode.mp_id_counter += 1
        self._id = "mp{}".format(MultiPathGraphNode.mp_id_counter)

        self._seq_length = self._compute_seq_length()

        self._component_id = 0  # will be set to a component ID after components are defined

        self._reweighted_flag = False # should be initialized to False before each reconstruction round.

    def get_id(self):
        return self._id
    
    def get_multiPathObj(self):
        return self._multiPath

    def get_multipath_graph(self):
        return self._mpg

    def get_simple_path(self):
        return(self._multiPath.get_simple_path())

    def set_reweighted_flag(self, flag_setting):
        self._reweighted_flag = flag_setting
        return

    def get_reweighted_flag(self):
        return self._reweighted_flag
            
        
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
        if self.get_reweighted_flag() is True:
            raise RuntimeError("Error, cant set weight on mpgn when reweighted flag is already set")

        self._prev_weight = self._weight
        
        self._weight = weight
        self.set_reweighted_flag(True)

        logger.debug("changed mpgn weight from {} to {}".format(self._prev_weight, self._weight))
        #sys.stderr.write("changed mpgn weight from {} to {}".format(self._prev_weight, self._weight))
        
        return


    def get_prev_weight(self):
        return self._prev_weight
    
    
    def get_score_EXCLUDE_containments(self, use_prev_weight=False):
        weight = self._prev_weight if use_prev_weight else self._weight
        score = self._count * weight
        if NORMALIZE_SCORES_BY_SEQ_LENGTH:
            score = score / self._seq_length  # normalize read count by feature length

        return score
    
    def get_score_INCLUDE_containments(self, use_prev_weight=False, mpgn_ignore=set()):
        
        seen = set(mpgn_ignore)

        total_counts = 0

        all_relevant_nodes = [self]
        contained_nodes = self.get_containments()
        if contained_nodes:
            all_relevant_nodes.extend(contained_nodes)

        for node in all_relevant_nodes:
            assert(type(node) == MultiPathGraphNode)
            if node not in seen:
                weight = node._prev_weight if use_prev_weight else node._weight
                total_counts += node._count * weight
        
        score = total_counts

        if NORMALIZE_SCORES_BY_SEQ_LENGTH:
            score = score / self._seq_length  # normalize read count by feature length
                
        return score
        
    
    def get_seq_length(self):
        return self._seq_length

    def get_component_id(self):
        return self._component_id

    def set_component_id(self, component_id):
        self._component_id = component_id
    
            
    def __repr__(self):
        return(self.toString(recursive=True)) # to get first round of containments


    def toString(self, recursive=False):
        containments = self.get_containments()
        
        text = "<mp:{} {}-{} Count:{} W:{:0.8f} Containments:{}, ScoreExcCont:{:.4f} ScoreInclCon:{:.4f} len:{}>".format(
            self.get_simple_path(),
            self._lend, self._rend, self._count, self._weight,
            len(containments),
            self.get_score_EXCLUDE_containments(use_prev_weight=False),
            self.get_score_INCLUDE_containments(use_prev_weight=False),
            self._seq_length)

        if recursive:
            for containment in containments:
                text += "\n\tcontained: " + containment.toString(recursive=False)

        return text

        
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


    def has_containments(self):
        if len(self._containments) > 0:
            return True
        else:
             return False
    
        
    def get_containments(self):
        return(list(self._containments))

    
    def contains_other_node(self, other_node):
                
        # this is spacer-aware
        if Simple_path_utils.simple_path_A_contains_and_compatible_with_simple_path_B(self.get_splice_graph(), self.get_simple_path(), other_node.get_simple_path()):
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
            if sg_node is not None:
                lend, rend = sg_node.get_coords()
                feature_len = rend - lend + 1
                if type(sg_node) == Exon:
                    seq_length += feature_len
                #elif type(sg_node) == Intron:
                #    seq_length += min(feature_len, 10000)  #//FIXME: determine max allowable length based on intron length distribution.
        
        return seq_length
        
    def reevaluate_weighting_via_path_compatibilities(self, transcript_multiPath):

        logger.debug("reevaluating weights for {}".format(self))
        #sys.stderr.write("reevaluating weights for {}".format(self))
        
        assert(type(transcript_multiPath) == MultiPath.MultiPath)
        
        transcript_simple_path = transcript_multiPath.get_simple_path()
        
        compatible_score = self.get_score_INCLUDE_containments(use_prev_weight=True)
        incompatible_score = 0

        connected_mpgns = self.get_predecessors() + self.get_successors()

        for node in connected_mpgns:
            node_simple_path = node.get_simple_path()
            if Simple_path_utils.path_A_contains_path_B(transcript_simple_path, node_simple_path):
                compatible_score += node.get_score_INCLUDE_containments(use_prev_weight=True)
            else:
                incompatible_score += node.get_score_INCLUDE_containments(use_prev_weight=True)

        pseudocount = 1e-3 
        
        fraction_compatible = (compatible_score + pseudocount) / (compatible_score + incompatible_score + pseudocount)
        
        current_weight = self.get_weight()

        adjusted_weight = current_weight - (fraction_compatible * current_weight)

        adjusted_weight = max(MIN_WEIGHT, adjusted_weight) # to avoid going too low.
        
        self.set_weight(adjusted_weight)

        logger.debug("fraction compatibile: {}, adjusted weight -> {}".format(fraction_compatible, adjusted_weight))
        
        return
        
        
