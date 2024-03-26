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

        def recursively_capture_nodes(mpgn):
            assert(type(mpgn) == MultiPathGraphNode)
            if mpgn not in self._all_represented_mpgns:
                self._all_represented_mpgns.add(mpgn)
                for contained_mpgn in mpgn.get_containments():
                    recursively_capture_nodes(contained_mpgn)
        
        for mpgn in path_list_of_multipath_graph_nodes:
            recursively_capture_nodes(mpgn)

                
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

        self._validate_compatible_containments()
        
        
        
    def __repr__(self):
        txt = "PASA_scored_path: {} (Score={:.5f}, InitScore={:.5f})\nmpgns:\n".format(self.get_multiPath_obj(), self.get_score(), self.get_initial_score())
        
        for mpgn in self.get_path_mpgn_list():
            txt += str(mpgn) + "\n"

        return txt
        
        
    def get_score(self):
        return self._score

    def get_initial_score(self):
        return self._initial_score

    def get_path_mpgn_list(self):
        return list(self._mpgn_list_path)

    def get_all_represented_mpgns(self, additional_mpgns_to_check=None):
        
        represented_mpgns = set(self._all_represented_mpgns)

        if additional_mpgns_to_check:

            scored_simple_path = self.get_multiPath_obj().get_simple_path()
            
            for mpgn in additional_mpgns_to_check:
                sg = mpgn.get_splice_graph()
                mpgn_simple_path = mpgn.get_simple_path()
                if Simple_path_utils.simple_path_A_contains_and_compatible_with_simple_path_B_spacer_aware_both_paths(sg, scored_simple_path, mpgn_simple_path):
                    represented_mpgns.add(mpgn)

        return represented_mpgns
    
    
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

        if self._score > self._initial_score:
            raise RuntimeError("Error, rescored path exceeds initial score for path: " + str(self))
        
            
        return
    

    def get_all_represented_reads(self):

        read_names = set()
        
        for mpgn in self._all_represented_mpgns:
            for read_name in mpgn.get_read_names():
                read_names.add(read_name)
            if mpgn.has_containments:
                for contained_mpgn in mpgn.get_containments():
                    for read_name in contained_mpgn.get_read_names():
                        read_names.add(read_name)
                        
        return read_names

    
    def toTranscript(self):

        mpgn_list = self.get_path_mpgn_list()

        mpgn_list = sorted(mpgn_list, key=lambda x: x._rend)
        
        splice_graph = mpgn_list[0].get_splice_graph()

        read_names = self.get_all_represented_reads()
        
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
        
        orient = '.' # '?'
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

        transcript_obj.add_read_names(read_names)
        
        return transcript_obj

    
        
    def compute_path_score(self):

        assert(self._cdna_len > 0 and self._contig_span_len > 0)
        
        score = 0

        mpgn_list = self.get_all_represented_mpgns()
        
        # sum up weights
        for mpgn in mpgn_list:
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
    
    

    def _validate_compatible_containments(self):

        
        scored_simple_path = self.get_multiPath_obj().get_simple_path()
        
        for mpgn in self.get_path_mpgn_list():
            sg = mpgn.get_splice_graph()
            mpgn_simple_path = mpgn.get_simple_path()
            
            if Simple_path_utils.simple_path_A_contains_and_compatible_with_simple_path_B_spacer_aware_both_paths(sg, scored_simple_path, mpgn_simple_path):
                logger.debug("Validated scored path: {}\ncontains path {}\n".format(scored_simple_path, mpgn_simple_path))
            else:
                raise RuntimeError("Error, scored path: {}\ndoes not contain path {}\n".format(scored_simple_path, mpgn_simple_path))

    
