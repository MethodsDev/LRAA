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
import PASA_SALRAA_Globals
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
from collections import defaultdict

logger = logging.getLogger(__name__)


class PASA_SALRAA:

    min_transcript_length = 200
    min_mpgn_read_count = 1

    max_contained_to_be_a_pasa_vertex = 10
    
    def __init__(self, splice_graph, num_parallel_processes=1):

        self._splice_graph = splice_graph

        self._multipath_graph = None  # set under build_multipath_graph()
        self._contig_acc = None # set under build_multipath_graph()
        
        self._num_parallel_processes = num_parallel_processes
        
        return


    def build_multipath_graph(self, contig_acc, contig_seq, bam_file, allow_spacers=False):

        logger.info(f"-building multipath graph for {contig_acc}")
        start_time = time.time()
        mp_counter = self._populate_read_multi_paths(contig_acc, contig_seq, bam_file)

        multipath_graph = MultiPathGraph(mp_counter, self._splice_graph, contig_acc, PASA_SALRAA.min_mpgn_read_count, allow_spacers)
        self._multipath_graph = multipath_graph
        self._contig_acc = contig_acc
        
        if PASA_SALRAA_Globals.DEBUG:
            ## debugging info
            logger.info("writing __multipath_graph.dat")
            multipath_graph.describe_graph("__multipath_graph.dat")
        
            
        build_time = time.time() - start_time
        logger.info("-multipath graph building took {:.1f} seconds.".format(build_time))
        
        return

    

    def reconstruct_isoforms(self, single_best_only=False):

        # define disjoint graph components.
        mpg = self._multipath_graph

        mpg_components = mpg.define_disjoint_graph_components_via_shared_splice_graph_vertex()

        logger.info("{} connected components identified".format(len(mpg_components)))
        mpg.write_mp_graph_nodes_to_gtf("__mpgns.pre.gtf")

        mpg_components = mpg.remove_small_components(mpg_components, PASA_SALRAA.min_transcript_length)
        logger.info("{} components surviving the min length {} criterion.".format(len(mpg_components), PASA_SALRAA.min_transcript_length))
        mpg.write_mp_graph_nodes_to_gtf("__mpgns.post_length_filter.gtf")
        
        all_reconstructed_transcripts = list()

        q = Queue()
        
        mpm = MultiProcessManager(self._num_parallel_processes, q)

        def get_mpgn_list_coord_span(mpgn_list):

            coords = list()
            for mpgn in mpgn_list:
                mpgn_coords = mpgn.get_coords()
                coords.extend(mpgn_coords)

            coords = sorted(coords)
            
            return coords[0], coords[-1]
        

        mpg_component_debug_dir = "__mpg_components"
        if PASA_SALRAA_Globals.DEBUG:
            if not os.path.exists(mpg_component_debug_dir):
                os.makedirs(mpg_component_debug_dir)

        def write_mpg_component_debug_file(mpgn_list, filename):
            logger.debug("-writing mpgn description file: {}".format(filename))
            with open(filename, 'wt') as ofh:
                for mpgn in mpgn_list:
                    print(str(mpgn), file=ofh)
        
        
        component_counter = 0
        for mpg_component in mpg_components:
            component_counter += 1
            coord_span = get_mpgn_list_coord_span(mpg_component)
            logger.info("PASA_SALRAA - assembly of component {} region: {}:{}-{}".format(component_counter, self._contig_acc, coord_span[0], coord_span[1]))
            mpg_token = "{}-{}-{}".format(self._contig_acc, coord_span[0], coord_span[1])
            if PASA_SALRAA_Globals.DEBUG:
                mpgn_description_filename = "{}/{}.mpgns.txt".format(mpg_component_debug_dir, mpg_token)
                write_mpg_component_debug_file(mpg_component, mpgn_description_filename)

                        
            p = Process(target=self._reconstruct_isoforms_single_component, name=mpg_token,
                        args=(q, mpg_component, component_counter, mpg_token, single_best_only) )

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
        
    
    def _reconstruct_isoforms_single_component(self, q, mpg_component, component_counter, mpg_token, single_best_only=False):


        ## exclude those mpgn's that are contained by many transcripts to reduce the trellis size.
        contained_mpg_counter = defaultdict(int)
        for mpgn in mpg_component:
            for mpgn_contained in mpgn.get_containments():
                contained_mpg_counter[mpgn_contained] += 1

        mpg_components_for_trellis = list()
        
        for mpgn in mpg_component:
            if contained_mpg_counter[mpgn] <= PASA_SALRAA.max_contained_to_be_a_pasa_vertex:
                mpg_components_for_trellis.append(mpgn)


        mpg_component = mpg_components_for_trellis # replace for trellis building
        logger.info("-num vertices for trellis: {}".format(len(mpg_component)))

        MIN_SCORE = PASA_SALRAA_Globals.config['min_path_score']

        best_transcript_paths = list()

        paths_seen = set()

        def reinit_weights(mpgn_list):
            for mpgn in mpgn_list:
                mpgn.set_reweighted_flag(False)
            return


        round_iter = 0        

        mpgns_require_representation = set(mpg_component)

        ###################
        ## Building trellis
        pasa_vertices = self._build_trellis(mpg_component, mpg_token)
        ###################

        
        if logger.getEffectiveLevel() == logging.DEBUG: ## for debugging info only
            if round_iter == 1:
                for pasa_vertex in pasa_vertices:
                    logger.debug(pasa_vertex.describe_pasa_vertex())

            self._write_all_scored_paths_to_file(component_counter, round_iter, mpg_token, pasa_vertices)

        all_scored_paths = self._retrieve_all_scored_paths(pasa_vertices)

        all_scored_paths = sorted(all_scored_paths, key=lambda x: x.get_score())
        
        while len(mpgns_require_representation) > 0 and len(all_scored_paths) > 0:

            round_iter += 1
            
            #reinit_weights(mpg_component)
            top_scored_path = all_scored_paths.pop()
            assert(type(top_scored_path) == PASA_scored_path)

            if top_scored_path.get_score() < MIN_SCORE:
                break
            
            mpgns_represented = top_scored_path.get_all_represented_mpgns(additional_mpgns_to_check=mpgns_require_representation)
            
            found_prev_unrepresented_mpgn = False
            for mpgn_represented in mpgns_represented:
                if mpgn_represented in mpgns_require_representation:
                    found_prev_unrepresented_mpgn = True
                    mpgns_require_representation.remove(mpgn_represented)

            if found_prev_unrepresented_mpgn:
                logger.debug("Retrieved best (Score={})  transcript path for mpg {} : {}".format(top_scored_path.get_score(), mpg_token, top_scored_path))
                self._write_best_score_path_info_to_file(top_scored_path, round_iter, mpg_token)
            
                best_transcript_paths.append(top_scored_path)

                # adjust weights
                self._decrement_transcript_path_vertices(top_scored_path, pasa_vertices)
                for path in all_scored_paths:
                    path.rescore()

                # reprioritize
                all_scored_paths = sorted(all_scored_paths, key=lambda x: x.get_score())
                

                
        # from the best transcript paths, reconstruct the actual transcripts themselves:

        self._validate_pairwise_incompatibilities(best_transcript_paths)
        
        
        transcripts = list() 


        for transcript_path in best_transcript_paths:
            assert(type(transcript_path) == PASA_scored_path)

            transcript_obj = transcript_path.toTranscript()

            #print(transcript_obj)
            if transcript_obj is not None:
                transcripts.append(transcript_obj)
                logger.debug("-assembled: {}".format(str(transcript_obj)))

        q.put(transcripts)
        

        

    def _populate_read_multi_paths(self, contig_acc, contig_seq, bam_file):

        """
        Reads the alignments from the BAM and for each read traces it
        through a path in the graph.
        The path is stored as a multipath object with a count associated with the number of reads assigned to it.
        """
        
        
        bam_extractor = Bam_alignment_extractor(bam_file)
        pretty_alignments = bam_extractor.get_read_alignments(contig_acc,
                                                              region_lend=self._splice_graph._region_lend,
                                                              region_rend=self._splice_graph._region_rend,
                                                              pretty=True)


        # grouping read alignments according to read pairings (for illumina PE data):
        # group alignments:  grouped_alignments['read_name'] = list(read1_pretty_alignment, read2_pretty_alignment, ...)
        grouped_alignments = self._group_alignments_by_read_name(pretty_alignments)
        
        logger.info("-got {} pretty alignments grouped into {} alignment groups.".format(len(pretty_alignments), len(grouped_alignments)))

        # distill read alignments into unique multipaths (so if 10k alignments yield the same structure, there's one multipath with 10k count associated)
        mp_counter = MultiPathCounter()

        # capture the read->path assignments:
        if PASA_SALRAA_Globals.DEBUG:
            read_graph_mappings_ofh = open("__read_graph_mappings.dat", "wt")
            

        logger.info("-start: mapping read alignments to the graph")
        for read_name in grouped_alignments:
            #print("{}\t{}".format(read_name, len(grouped_alignments[read_name])))
            paths_list = list()
            read_type = None
            for pretty_alignment in grouped_alignments[read_name]:
                if read_type is None:
                    read_type = pretty_alignment.get_read_type() 

                path = self._map_read_to_graph(pretty_alignment.get_pretty_alignment_segments())
                                    
                logger.debug("pretty_alignment: {} maps to graph path: {}".format(pretty_alignment, path))
                if path and path != SPACER:
                    assert path[0] != SPACER, "path[0] is SPACER, not allowed"
                    assert path[-1] != SPACER, "path[-1] is SPACER, not allowed"
                    paths_list.append(path)

            mp = None
            if paths_list:
                mp = MultiPath(self._splice_graph, paths_list, read_types = { read_type, }, read_names = { read_name, } )
                logger.debug("paths_list: {} -> mp: {}".format(paths_list, mp))
                mp_counter.add(mp)

            if PASA_SALRAA_Globals.DEBUG:
                read_graph_mappings_ofh.write("\t".join([read_name, str(pretty_alignment), str(mp)]) + "\n")


        if PASA_SALRAA_Globals.DEBUG:
            read_graph_mappings_ofh.close()
        #print(mp_counter)

        logger.info("-done: mapping read alignments to the graph")
            
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
                if not terminal_segment and path_part != [SPACER]:
                    terminal_segment = [SPACER]
                if terminal_segment:
                    path_part.extend(terminal_segment)
            else:
                # internal segment
                #   first, get preceding intron
                path_part = self._get_intron_node_id(alignment_segments[i-1], segment)
                if not path_part:
                    path_part = [SPACER]
                internal_segment = self._map_segment_to_graph_INTERNAL(segment)
                if not internal_segment and path_part != [SPACER]:
                    internal_segment = [SPACER]
                if internal_segment:
                    path_part.extend(internal_segment)

            logger.debug("segment: {}  mapped to {}".format(segment, path_part))
                    
            if path_part:
                path.extend(path_part)
                #print("\tpath extended to: {}".format(path))
            else:
                if len(path) == 0 or path[-1] != SPACER:
                    path.append(SPACER) # spacer



        # refine path

        if SPACER in path:
            path = Simple_path_utils.trim_terminal_spacers(path)
            path = self._remove_stutters(path)

        path = Simple_path_utils.add_spacers_between_disconnected_nodes(self.get_splice_graph(), path)
        
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
                exon_segment._rend <= segment[1]):

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
                exon_segment._lend >= segment[0]):

                path.append(exon_segment.get_id())

        #logger.debug("segment {} maps to {}".format(segment, path))
                
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


    def _remove_stutters(self, path):

        new_path = list()
        
        for i, node_id in enumerate(path):
            if node_id != SPACER:
                if len(new_path) == 0:
                    new_path.append(node_id)
                elif new_path[-1] != SPACER:
                    if new_path[-1] != node_id:
                        new_path.append(node_id)
                elif new_path[-1] == SPACER and len(new_path) > 1 and new_path[-2] != node_id:
                    new_path.append(node_id)
            # it is a spacer
            elif len(new_path) > 0:
                if new_path[-1] != SPACER:
                    new_path.append(SPACER)

        new_path = Simple_path_utils.trim_terminal_spacers(new_path)

        if new_path != path:
            logger.debug("{} -> {}".format(path, new_path))
        
        return(new_path)
    
        

    def _try_easy_fill_spacers(self, path):

        return Simple_path_utils.try_fill_spacers_via_splicegraph(self._splice_graph, path)
        
    
    
    def _build_trellis(self, mpg_component, mpg_token):

        logger.debug("-building trellis for {}".format(mpg_token))
        
        mpg = self._multipath_graph

        nodes = sorted(mpg_component, key=lambda x: (x._lend, x._rend))
        
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

        logger.debug("-retrieved best transcript path {} with score: {}".format(best_scoring_path, best_score))
                    
        return best_scoring_path


    def _retrieve_all_scored_paths(self, pasa_vertices):

        all_scored_paths = list()
        
        for pasa_vertex in pasa_vertices:
            these_scored_paths = pasa_vertex.get_fromPaths()
            all_scored_paths.extend(these_scored_paths)
                    
        return all_scored_paths


    
    

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


        def recursive_reweight(mpgn):
            if mpgn.get_reweighted_flag() is False:
                mpgn.set_weight(0.000001)
                for contained_mpgn in mpgn.get_containments():
                    recursive_reweight(contained_mpgn)
        
        for mpgn in mpgn_list:
            logger.debug("_decrement: {}".format(mpgn))
            recursive_reweight(mpgn)
            
            #if mpgn.get_reweighted_flag() is False:
            #    mpgn.reevaluate_weighting_via_path_compatibilities(transcript_path_multipath_obj)

        ## //TODO: //FIXME  why vertex path incompatible but still part of transcript path?

        for mpgn in transcript_path.get_path_mpgn_list():
            recursive_reweight(mpgn)
        
        return


    def _rescore_transcript_paths(self, pasa_vertices):
        
        for pasa_vertex in pasa_vertices:
            pasa_vertex.rescore_fromPaths()

        return
    
    def _write_all_scored_paths_to_file(self, component_counter, round_iter, mpg_token, pasa_vertices):

        
        outdirname = "__all_scored_paths"

        if not os.path.exists(outdirname):
            os.makedirs(outdirname)

        outputfilename = "{}/scored_paths.{}.R{}.gtf".format(outdirname, mpg_token, round_iter)
        
        ofh = open(outputfilename, 'wt')

        for pasa_vertex in pasa_vertices:
            print ("## pasa vertex:", file=ofh)
            print(pasa_vertex.describe_pasa_vertex(), file=ofh)
            from_paths = pasa_vertex.get_fromPaths()
            for from_path in from_paths:
                trans_obj = from_path.toTranscript()
                trans_obj.add_meta('score', from_path.get_score())
                gtf = trans_obj.to_GTF_format()
                ofh.write(gtf + "\n")
            print("", file=ofh) # spacer"

        ofh.close()

        return
    

    def _write_best_score_path_info_to_file(self, transcript_path, round_iter, mpg_token):

        outdirname = "__selected_best_scored_paths"
        
        if not os.path.exists(outdirname):
            os.makedirs(outdirname)

        outputfilename = "{}/selected_best_path.{}.R{}.gtf".format(outdirname, mpg_token, round_iter)
        
        with open(outputfilename, 'wt') as ofh:
            print("Transcript path: " + str(transcript_path), file=ofh)

 
    
    def get_splice_graph(self):
        return self._splice_graph



    def _validate_pairwise_incompatibilities(self, pasa_scored_paths):

        if len(pasa_scored_paths) < 2:
            return
        
        for i in range(1, len(pasa_scored_paths)):
            mp_A = pasa_scored_paths[i].get_multiPath_obj()  
            sp_A = mp_A.get_simple_path()
            sg = mp_A.get_splice_graph()
            
            for j in range(i-1, -1, -1):
                mp_B = pasa_scored_paths[j].get_multiPath_obj()
                sp_B = mp_B.get_simple_path()

                if Simple_path_utils.simple_paths_overlap_and_compatible_spacer_aware_both_paths(sg, sp_A, sp_B):
                    raise RuntimeError("\n\n\n##****************************************************************************\n" +
                                       "Error, final paths:\n{} and \n{}\noverlap and are compatible - should have gotten assembled together".format(pasa_scored_paths[i], pasa_scored_paths[j]) + 
                                       "\n##***************************************************************************************\n\n")
                
        return
    
                
