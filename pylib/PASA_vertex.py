import sys, os, re
from MultiPath import MultiPath
from MultiPathGraph import MultiPathGraphNode
from PASA_scored_path import PASA_scored_path

import logging

logger = logging.getLogger(__name__)

class PASA_vertex:


    def __init__(self, multipath_graph_node):

        assert(type(multipath_graph_node) == MultiPathGraphNode)
        
        self._multipath_graph_node = multipath_graph_node            
        self._fromPaths = list() # hold scored paths.

        # add unextended current path node as initial path

        aggregate_counts = multipath_graph_node.get_count()
        for containment_node in multipath_graph_node.get_containments():
            aggregate_counts += containment_node.get_count()

        # score = reads normalized by seq length
        mpgn_seq_len = multipath_graph_node.get_seq_length()
        print("mpgn: {}".format(multipath_graph_node))
        init_score = aggregate_counts / mpgn_seq_len
        initial_scored_path = PASA_scored_path([self.get_mpgn()], init_score)
        
        self._fromPaths.append(initial_scored_path)


        return

    
    def __repr__(self):
        return("PASA_vertex for {}".format(self._multipath_graph_node))
        
    
    def get_multipath_graph_node(self):
        return self._multipath_graph_node

    def get_mpgn(self):
        return self.get_multipath_graph_node()
    

    def get_multipath_graph(self):
        return self._multipath_graph_node.get_multipath_graph()

    
    def get_fromPaths(self):
        return(list(self._fromPaths))

    
    def add_highest_scoring_path_extension(self, prev_pasa_vertex):

        best_score = 0
        best_prev_scored_path = None

        self_mpgn = self.get_mpgn()
        assert(type(self_mpgn) == MultiPathGraphNode)
        
        for prev_scored_path in prev_pasa_vertex.get_fromPaths():
            if not prev_scored_path.incompatibility_detected(self_mpgn):

                extension_path_candidate = prev_scored_path.create_scored_path_extension(self_mpgn)
                if extension_path_candidate.get_score() > best_score:
                    best_score = extension_path_candidate.get_score()
                    best_prev_scored_path = extension_path_candidate
                    
        if best_prev_scored_path is not None:
            self._fromPaths.append(best_prev_scored_path)
            logger.debug("Added best extension path: {} to {}".format(best_prev_scored_path, self))
            
        return
    

    def rescore_fromPaths(self):

        for scored_path in self.get_fromPaths():
            scored_path.rescore()

        return
        


    def describe_pasa_vertex(self):

        pv_text = str(self) + "\n"

        from_paths = self.get_fromPaths()

        from_paths = sorted(from_paths, key=lambda x: x._score, reverse=True)

        for path in from_paths:
            pv_text += "\t" + str(path) + "\n"

        return pv_text
    
