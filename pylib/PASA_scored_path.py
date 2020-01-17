import sys, os, re
from MultiPath import MultiPath
from MultiPathGraph import MultiPathGraphNode

import math
import logging

logger = logging.getLogger(__name__)

class PASA_scored_path:

    def __init__(self, path_list_of_multipath_graph_nodes):


        for mpgn in path_list_of_multipath_graph_nodes:
            assert(type(mpgn) == MultiPathGraphNode)
        
        self._mpgn_list_path = path_list_of_multipath_graph_nodes
        
        self._multiPath_obj = MultiPath.multiPath_from_mpgn_list(self._mpgn_list_path)

        self._cdna_len = self._multiPath_obj.get_cdna_length()

        span_lend, span_rend = self._multiPath_obj.get_coords()
        
        self._contig_span_len = span_rend - span_lend + 1
        
        score = self.compute_path_score()
        
        self._score = score
        self._initial_score = score

        
        
    def __repr__(self):
        return("PASA_scored_path: (score={}, IScore={}) mpgns: {}".format(self.get_score(), self.get_initial_score(), self.get_path_mpgn_list()))
    
        
    def get_score(self):
        return self._score

    def get_initial_score(self):
        return self._initial_score

    def get_path_mpgn_list(self):
        return list(self._mpgn_list_path)

    
    def get_multiPath_obj(self):
        return self._multiPath_obj
    
    
    def incompatibility_detected(self, extension_mpgn):

        mpg = extension_mpgn.get_multipath_graph()

        for mpgn in self.get_path_mpgn_list():
            if mpg.incompatible_mpgn_pair(mpgn, extension_mpgn):
                return True

        return False

    def create_scored_path_extension(self, extension_mpgn):

        path_list = self.get_path_mpgn_list() + [extension_mpgn]

        extension_scored_path = PASA_scored_path(path_list)

        return extension_scored_path


    def rescore(self):
        self._score = self.compute_path_score()
        return
    
    
    def toTranscript(self):

        mpgn_list = self.get_path_mpgn_list()

        mpgn_list = sorted(mpgn_list, key=lambda x: x._rend)
        
        splice_graph = mpgn_list[0].get_splice_graph()
        
        # merge to a single multipath object
        simple_path_list = list()
        for mpgn in mpgn_list:
            mp = mpgn.get_multiPathObj()
            simple_path_list.append(mp.get_simple_path())
        
        transcript_mp = MultiPath(splice_graph, simple_path_list)

        return transcript_mp

    
        
    def compute_path_score(self):

        assert(self._cdna_len > 0 and self._contig_span_len > 0)
        
        score = 0
        seen = set()


        # ordered list, containments only incorporated into score accoridng to left-most containment node.
        mpgn_list = self.get_path_mpgn_list()
        mpgn_list = sorted(mpgn_list, key=lambda x: x._lend)

        for mpgn in mpgn_list:
            if mpgn not in seen:
                score += mpgn.get_score_INCLUDE_containments(use_prev_weight=False, mpgn_ignore=seen)
            
            seen.add(mpgn)
            for containment in mpgn.get_containments():
                if containment not in seen:
                    seen.add(containment)
        
        return score
    
    
