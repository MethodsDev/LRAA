import sys, os, re
from MultiPath import MultiPath
from MultiPathGraph import MultiPathGraphNode
from GenomeFeature import *
import Simple_path_utils
from Transcript import Transcript
import PASA_SALRAA_Globals

import math
import logging

logger = logging.getLogger(__name__)

class PASA_scored_path:

    def __init__(self, path_list_of_multipath_graph_nodes):


        self._all_represented_mpgns = set() # stores all input mpgns and their contained mpgns
        
        for mpgn in path_list_of_multipath_graph_nodes:
            assert(type(mpgn) == MultiPathGraphNode)
            self._all_represented_mpgns.add(mpgn)
            for mpgn_contained in mpgn.get_containments():
                self._all_represented_mpgns.add(mpgn_contained)
        
        self._mpgn_list_path = path_list_of_multipath_graph_nodes
        
        self._multiPath_obj = MultiPath.multiPath_from_mpgn_list(self._mpgn_list_path)

        self._cdna_len = self._multiPath_obj.get_cdna_length()

        span_lend, span_rend = self._multiPath_obj.get_coords()
        
        self._contig_span_len = span_rend - span_lend + 1

        # init
        self._score = -1
        self._initial_score = -1
        
        score = self.compute_path_score()
        
        # set
        self._score = score
        self._initial_score = score

        
        
    def __repr__(self):
        return("PASA_scored_path: (score={:.5f}, IScore={:.5f}) mpgns: {}".format(self.get_score(), self.get_initial_score(), self.get_path_mpgn_list()))


        
        
    def get_score(self):
        return self._score

    def get_initial_score(self):
        return self._initial_score

    def get_path_mpgn_list(self):
        return list(self._mpgn_list_path)

    def get_all_represented_mpgns(self):
        return list(self._all_represented_mpgns)
    
    
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

        exons_and_introns = transcript_mp.get_ordered_exons_and_introns()

        #print("exons and introns: ")
        #for exon_or_intron in exons_and_introns:
        #    print("\t" + str(exon_or_intron))
        
        transcript_exon_segments = list()

        orient = '?'
        contig_acc = exons_and_introns[0].get_contig_acc()

        for feature in exons_and_introns:
            if type(feature) == Exon:
                transcript_exon_segments.append(feature.get_coords())
            elif type(feature) == Intron:
                orient = feature.get_orient()

        if len(transcript_exon_segments) == 0:
            logger.warning("bug - shouldn't have exonless transcript features: {}".format(transcript_path)) # //FIXME: bug
            return None

        transcript_exon_segments = Simple_path_utils.merge_adjacent_segments(transcript_exon_segments)

        #print("merged segments: " + str(transcript_exon_segments))

        transcript_obj = Transcript(contig_acc, transcript_exon_segments, orient)

        #transcript_obj.set_scored_path_obj(self)  # commenting out / lightening up object
        
        return transcript_obj

    
        
    def compute_path_score(self):

        assert(self._cdna_len > 0 and self._contig_span_len > 0)
        
        score = 0


        all_mpgn_nodes = set()
        
        
        mpgn_list = self.get_path_mpgn_list()


        def recursively_capture_nodes(mpgn):

            if mpgn not in all_mpgn_nodes:
                all_mpgn_nodes.add(mpgn)
                for contained_mpgn in mpgn.get_containments():
                    recursively_capture_nodes(mpgn)

        for mpgn in mpgn_list:
            recursively_capture_nodes(mpgn)

        # sum up weights
        for mpgn in all_mpgn_nodes:
            score += mpgn.get_weight() * mpgn.get_count()

        """


        audit_txt = ""
        
        if PASA_SALRAA_Globals.DEBUG:
            audit_txt = "Computing path score for: {}\n".format(self._multiPath_obj)
        
        for mpgn in mpgn_list:
            if mpgn not in seen:
                score_increment = mpgn.get_score_INCLUDE_containments(use_prev_weight=False, mpgn_ignore=seen)
                score += score_increment

                if PASA_SALRAA_Globals.DEBUG:
                    audit_txt += "\tscore subtotal: {:.5f}, increment {:.5f}, mpgn: {}\n".format(score, score_increment, mpgn)
                
            seen.add(mpgn)
            for containment in mpgn.get_containments():
                if containment not in seen:
                    seen.add(containment)
                    if PASA_SALRAA_Globals.DEBUG:
                        audit_txt += "\t\tcontainment: {}\n".format(containment)


        if PASA_SALRAA_Globals.DEBUG:
            logger.debug(audit_txt)
           
        """

        logger.debug(str(self) + "\n^^^ computed with score = {}".format(score))
            
        return score
    
    
