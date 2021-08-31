#!/usr/bin/env python
# encoding: utf-8

import sys, os, re
import subprocess
import logging
import string
import pysam
from collections import defaultdict
import networkx as nx
import intervaltree as itree
from GenomeFeature import *
from Bam_alignment_extractor import Bam_alignment_extractor
from MultiPath import MultiPath
from MultiPathCounter import MultiPathCounter
from PASA_SALRAA_Globals import SPACER
from MultiPathGraph import MultiPathGraph
from PASA_vertex import PASA_vertex
from Transcript import Transcript
import Simple_path_utils
from PASA_scored_path import PASA_scored_path
from shutil import rmtree
import time
from multiprocessing import Process, Queue
import traceback
from MultiProcessManager import MultiProcessManager

logger = logging.getLogger(__name__)


DEBUG_FILES_FLAG = False

class PASA_SALRAA:

    min_transcript_length = 200
    min_mpgn_read_count = 1

    def __init__(self, splice_graph, num_parallel_processes=1):

        self._splice_graph = splice_graph

        self._multipath_graph = None  # set under build_multipath_graph()

        self._num_parallel_processes = num_parallel_processes
        
        return


    def build_multipath_graph(self, contig_acc, contig_seq, bam_file):

        mp_counter = self._populate_read_multi_paths(contig_acc, contig_seq, bam_file)

        multipath_graph = MultiPathGraph(mp_counter, self._splice_graph, PASA_SALRAA.min_mpgn_read_count)
        self._multipath_graph = multipath_graph

        ## debugging info
        logger.info("writing __multipath_graph.dat")
        multipath_graph.describe_graph("__multipath_graph.dat")

        
        return

    

    def reconstruct_isoforms(self, single_best_only=False):

        # define disjoint graph components.
        mpg = self._multipath_graph

        mpg_components = mpg.define_disjoint_graph_components_via_shared_splice_graph_vertex()

        logger.info("{} connected components identified".format(len(mpg_components)))
        mpg.write_mp_graph_nodes_to_gtf("mpgns.pre.gtf")

        mpg_components = mpg.remove_small_components(mpg_components, PASA_SALRAA.min_transcript_length)
        logger.info("{} components surviving the min length {} criterion.".format(len(mpg_components), PASA_SALRAA.min_transcript_length))
        mpg.write_mp_graph_nodes_to_gtf("mpgns.post_length_filter.gtf")
        
        all_reconstructed_transcripts = list()

        q = Queue()
        
        mpm = MultiProcessManager(self._num_parallel_processes, q)
        
        component_counter = 0
        for mpg_component in mpg_components:
            component_counter += 1
            logger.info("PASA_SALRAA - assembly of component {}".format(component_counter))
            
            p = Process(target=self._reconstruct_isoforms_single_component,
                        args=(q, mpg_component, component_counter, single_best_only) )

            mpm.launch_process(p)

        num_failures = mpm.wait_for_remaining_processes()

        logger.info(mpm.summarize_status())
        
        if num_failures:
            raise RuntimeError("Error, {} component failures encountered".format(num_failures))

        all_reconstructed_transcripts = list()

        queue_contents = mpm.retrieve_queue_contents()
        for entry in queue_contents:
            all_reconstructed_transcripts.extend(entry)
                            
                
        return all_reconstructed_transcripts


    
    ##################
    ## Private methods
    ##################
        
    
    def _reconstruct_isoforms_single_component(self, q, mpg_component, component_counter, single_best_only=False):


        MIN_SCORE_RATIO = 0.0001

        best_transcript_paths = list()

        paths_seen = set()


        def reinit_weights(mpgn_list):
            for mpgn in mpgn_list:
                mpgn.set_reweighted_flag(False)
            return


        round_iter = 0        
        while True:

            round_iter += 1

            reinit_weights(mpg_component)

            pasa_vertices = self._build_trellis(mpg_component)

            if logger.getEffectiveLevel() == logging.DEBUG: ## for debugging info only
                if round_iter == 1:
                    for pasa_vertex in pasa_vertices:
                        logger.debug(pasa_vertex.describe_pasa_vertex())

                self._write_all_scored_paths_to_file(component_counter, round_iter, pasa_vertices)



            transcript_path = self._retrieve_best_transcript(pasa_vertices)

            if transcript_path is None:
                break

            assert(type(transcript_path) == PASA_scored_path)

            logger.debug("Retrieved best transcript path: {}".format(transcript_path))

            transcript_path_token = str(transcript_path.get_multiPath_obj())
            if transcript_path_token in paths_seen:
                logger.debug("best path {} already reported. Stopping path extractions from component now.".format(transcript_path_token))
                break

            paths_seen.add(transcript_path_token)


            if (transcript_path.get_score() > 0 and
                (len(best_transcript_paths) == 0 or
                 transcript_path.get_score() / best_transcript_paths[0].get_initial_score() >= MIN_SCORE_RATIO) ):

                best_transcript_paths.append(transcript_path)

                if single_best_only:
                    break

                self._decrement_transcript_path_vertices(transcript_path, pasa_vertices)
                self._rescore_transcript_paths(pasa_vertices)
            else:
                break


        # from the best transcript paths, reconstruct the actual transcripts themselves:


        transcripts = list() 


        for transcript_path in best_transcript_paths:
            assert(type(transcript_path) == PASA_scored_path)

            transcript_obj = transcript_path.toTranscript()

            #print(transcript_obj)
            if transcript_obj is not None:
                transcripts.append(transcript_obj)

        q.put(transcripts)
        

        

    def _populate_read_multi_paths(self, contig_acc, contig_seq, bam_file):

        bam_extractor = Bam_alignment_extractor(bam_file)
        pretty_alignments = bam_extractor.get_read_alignments(contig_acc, pretty=True)

        grouped_alignments = self._group_alignments_by_read_name(pretty_alignments)

        logger.info("-got {} pretty alignments grouped into {} alignment groups.".format(len(pretty_alignments), len(grouped_alignments)))
        
        mp_counter = MultiPathCounter()

        read_graph_mappings_ofh = open("__read_graph_mappings.dat", "wt")
        
        for read_name in grouped_alignments:
            #print("{}\t{}".format(read_name, len(grouped_alignments[read_name])))
            paths_list = list()
            for pretty_alignment in grouped_alignments[read_name]:
                path = self._map_read_to_graph(pretty_alignment.get_pretty_alignment_segments())
                #print("pretty_alignment: {} maps to graph path: {}".format(pretty_alignment, path))
                if path and path != SPACER:
                    paths_list.append(path)

            mp = None
            if paths_list:
                mp = MultiPath(self._splice_graph, paths_list)
                #print("paths_list: {} -> mp: {}".format(paths_list, mp))
                mp_counter.add(mp)

            read_graph_mappings_ofh.write("\t".join([read_name, str(pretty_alignment), str(mp)]) + "\n")


        read_graph_mappings_ofh.close()
        #print(mp_counter)
        
        return mp_counter
    
         

    def _group_alignments_by_read_name(self, pretty_alignments):

        grouped_alignments = defaultdict(list)

        for pretty_alignment in pretty_alignments:
            pysam_alignment = pretty_alignment.get_pysam_alignment()
            read_name = pysam_alignment.query_name
            grouped_alignments[read_name].append(pretty_alignment)

        return grouped_alignments
    

    def _map_read_to_graph(self, alignment_segments):

        path = list()

        num_segments = len(alignment_segments)

        for i in range(num_segments):

            segment = alignment_segments[i]

            path_part = None
            
            ## determine type of segment
            if i == 0 and num_segments == 1:
                # single exon segment type
                path_part = self._map_segment_to_graph_SINGLE(segment)
            elif i == 0:
                # initial segment
                path_part = self._map_segment_to_graph_INITIAL(segment)
            elif i == num_segments - 1:
                # terminal segment
                path_part = self._get_intron_node_id(alignment_segments[i-1], segment)
                if not path_part:
                    path_part = [SPACER]
                terminal_segment = self._map_segment_to_graph_TERMINAL(segment)
                if not terminal_segment:
                    terminal_segment = [SPACER]
                path_part.extend(terminal_segment)
            else:
                # internal segment
                #   first, get preceding intron
                path_part = self._get_intron_node_id(alignment_segments[i-1], segment)
                if not path_part:
                    path_part = [SPACER]
                internal_segment = self._map_segment_to_graph_INTERNAL(segment)
                if not internal_segment:
                    internal_segment = [SPACER]
                path_part.extend(internal_segment)

            print("segment: {}  mapped to {}".format(segment, path_part))
                    
            if path_part:
                path.extend(path_part)
                #print("\tpath extended to: {}".format(path))
            else:
                if len(path) == 0 or path[-1] != SPACER:
                    path.append(SPACER) # spacer


        # trim any terminal spacer
        while len(path) > 0  and path[-1] == SPACER:
            path = path[:-1]

        # trim any initial spacer
        while len(path) > 0 and path[0] == SPACER:
            path.pop(0)
            

        if SPACER in path:
            path = self._try_easy_fill_spacers(path)
            
            
        return path
    

    def _get_intron_node_id(self, prev_segment, next_segment):

        intron_lend = prev_segment[1] + 1
        intron_rend = next_segment[0] - 1

        intron_obj = self._splice_graph.get_intron_node_obj(intron_lend, intron_rend)
        if intron_obj:
            return [intron_obj.get_id()]
        else:
            return None


    def _map_segment_to_graph_SINGLE(self, segment):

        overlapping_segments = self._splice_graph.get_overlapping_exon_segments(segment[0], segment[1])

        overlapping_segments = sorted(overlapping_segments, key=lambda x: x._lend)

        path = list()
        for exon_segment in overlapping_segments:
            id = exon_segment.get_id()
            path.append(id)

        return path

    def _map_segment_to_graph_INITIAL(self, segment):

        path = list()

        overlapping_segments = self._splice_graph.get_overlapping_exon_segments(segment[0], segment[1])
        overlapping_segments = sorted(overlapping_segments, key=lambda x: x._lend)
        
        for exon_segment in overlapping_segments:
            # check for overlap and not extending beyond feature rend
            if (segment[0] < exon_segment._rend and
                segment[1] > exon_segment._lend and
                segment[1] <= exon_segment._rend):

                path.append(exon_segment.get_id())

        return path


    
    def _map_segment_to_graph_TERMINAL(self, segment):

        path = list()

        overlapping_segments = self._splice_graph.get_overlapping_exon_segments(segment[0], segment[1])
        overlapping_segments = sorted(overlapping_segments, key=lambda x: x._lend)
        
        for exon_segment in overlapping_segments:
            # check for overlap and not extending beyond feature rend
            if (segment[0] < exon_segment._rend and
                segment[1] > exon_segment._lend and
                segment[0] >= exon_segment._lend):

                path.append(exon_segment.get_id())

        return path


    def _map_segment_to_graph_INTERNAL(self, segment):
        
        path = list()

        overlapping_segments = self._splice_graph.get_overlapping_exon_segments(segment[0], segment[1])
        overlapping_segments = sorted(overlapping_segments, key=lambda x: x._lend)
        
        for exon_segment in overlapping_segments:
            # check for overlap only  (prev: and not extending beyond feature rend)
            if (segment[0] < exon_segment._rend and
                segment[1] > exon_segment._lend):
                #and
                #segment[0] <= exon_segment._lend and
                #exon_segment._rend <= segment[1]):

                path.append(exon_segment.get_id())

        return path


    def _try_easy_fill_spacers(self, path):


        new_path = path.copy()
        
        for i in range(len(path)):
            if path[i] == SPACER:
                if i == 0 or i == len(path)-1:
                    # cant be terminal
                    continue
                
                prev_node = path[i-1]
                next_node = path[i+1]

                if prev_node == SPACER or next_node == SPACER:
                    continue
                
                prev_node_obj = self._splice_graph.get_node_obj_via_id(prev_node)
                next_node_obj = self._splice_graph.get_node_obj_via_id(next_node)

                # require exons on both sides of spacer
                if not ( type(prev_node_obj) == Exon and type(next_node_obj) == Exon):
                    continue
                
                intron_lend = prev_node_obj._rend + 1
                intron_rend = next_node_obj._lend - 1

                intron_obj = self._splice_graph.get_intron_node_obj(intron_lend, intron_rend)
                if intron_obj:
                    intron_id = intron_obj.get_id()
                    new_path[i] = intron_id # replace spacer

        
        return new_path
                

    
    
    def _build_trellis(self, mpg_component):

        mpg = self._multipath_graph

        nodes = sorted(mpg_component, key=lambda x: x._lend)
        
        # init the pasa vertex list
        
        pasa_vertices = list()
        
        for node in nodes:
            pasa_vertex = PASA_vertex(node)
            pasa_vertices.append(pasa_vertex)

        
        for i in range(1, len(pasa_vertices)):
            pv_i = pasa_vertices[i]

            for j in range(i-1, -1, -1):
                pv_j = pasa_vertices[j]
                
                if mpg.has_edge(pv_j.get_multipath_graph_node(), pv_i.get_multipath_graph_node()):
                    pv_i.add_highest_scoring_path_extension(pv_j)
                    

        return pasa_vertices


    def _retrieve_best_transcript(self, pasa_vertices):

        best_scoring_path = None
        best_score = 0
        
        for pasa_vertex in pasa_vertices:
            
            scored_paths = pasa_vertex.get_fromPaths()
            for scored_path in scored_paths:
                if scored_path.get_score() > best_score:
                    best_score = scored_path.get_score()
                    best_scoring_path = scored_path

        return best_scoring_path


    def _decrement_transcript_path_vertices(self, transcript_path, pasa_vertices):

        logger.debug("_decrement_transcript_path_vertices")
        
        assert(type(transcript_path) == PASA_scored_path)

        # examine all pasa vertices that are contained and compatible with the transcript_path 
        mpgn_list = list()
        
        transcript_path_multipath_obj = transcript_path.get_multiPath_obj()

        mpgns_not_compatible = list()
        
        for pasa_vertex in pasa_vertices:
            mpgn = pasa_vertex.get_multipath_graph_node()
            other_multipath_obj = mpgn.get_multiPathObj()
            if transcript_path_multipath_obj.is_overlapping_contained_and_compatible(other_multipath_obj):
                mpgn_list.append(mpgn)
            else:
                mpgns_not_compatible.append(mpgn)
                
        logger.debug("mpgns found compatible with transcript path: {} include {}".format(transcript_path_multipath_obj, mpgn_list))
        logger.debug("mpgns found INcomptable are: {}".format(mpgns_not_compatible))
        
        for mpgn in mpgn_list:
            logger.debug("_decrement: {}".format(mpgn))
            if mpgn.get_reweighted_flag() is False:
                mpgn.reevaluate_weighting_via_path_compatibilities(transcript_path_multipath_obj)
                
                        
        return


    def _rescore_transcript_paths(self, pasa_vertices):
        
        for pasa_vertex in pasa_vertices:
            pasa_vertex.rescore_fromPaths()

        return
    
    def _write_all_scored_paths_to_file(self, component_counter, round_iter, pasa_vertices):

        global DEBUG_FILES_FLAG
        
        outdirname = "__all_scored_paths"
        if DEBUG_FILES_FLAG is False:
            if os.path.exists(outdirname):
                rmtree(outdirname)
                
            os.makedirs(outdirname)
            DEBUG_FILES_FLAG = True

        outputfilename = "{}/scored_paths.R{}.gtf".format(outdirname, round_iter)
        
        mode = 'at' if os.path.exists(outputfilename) else 'wt'
        ofh = open(outputfilename, mode)

        for pasa_vertex in pasa_vertices:
            from_paths = pasa_vertex.get_fromPaths()
            for from_path in from_paths:
                trans_obj = from_path.toTranscript()
                trans_obj.add_meta('score', from_path.get_score())
                gtf = trans_obj.to_GTF_format()
                ofh.write(gtf + "\n")

        ofh.close()

        return
    
        
def ladeda(q, stuff):
    q.put("yowzer!")
    
