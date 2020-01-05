import sys, os, re
from MultiPathGraph import MultiPathGraphNode


class PASA_scored_path:

    def __init__(self, path_list_of_multipath_graph_nodes, score):


        for mpgn in path_list_of_multipath_graph_nodes:
            assert(type(mpgn) == MultiPathGraphNode)
        
        self._path = path_list_of_multipath_graph_nodes
        self._score = score
        

    def get_score(self):
        return self._score

    
    def incompatibility_detected(self, extension_mpgn):

        mpg = extension_mpgn.get_multipath_graph()

        for mpgn in self._path:
            if mpg.incompatible_mpgn_pair(mpgn, extension_mpgn):
                return True

        return False

    def create_scored_path_extension(self, extension_mpgn):

        score = 0

        path_list = self._path + [extension_mpgn]

        seen = set()
        for mpgn in path_list:
            score = mpgn.get_count()
            seen.add(mpgn)
            for containment in mpgn.get_containments():
                if containment not in seen:
                    seen.add(containment)
                    score += containment.get_count()

        extension_scored_path = PASA_scored_path(path_list, score)

        return extension_scored_path

    

class PASA_vertex:


    def __init__(self, multipath_graph_node):

        assert(type(multipath_graph_node) == MultiPathGraphNode)
        
        self._multipath_graph_node = multipath_graph_node            
        self._fromPaths = list() # hold scored paths.

        # add unextended current path node as initial path

        init_score = multipath_graph_node.get_count()
        for containment_node in multipath_graph_node.get_containments():
            init_score += containment_node.get_count()
        
        initial_scored_path = PASA_scored_path([self.get_mpgn()], init_score)

        self._fromPaths.append(initial_scored_path)


        return


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
                
        return
    
