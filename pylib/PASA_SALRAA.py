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

logger = logging.getLogger(__name__)

class PASA_SALRAA:


    def __init__(self, splice_graph):

        self._splice_graph = splice_graph

        return


    def populate_read_multi_paths(self, contig_acc, contig_seq, bam_file):

        bam_extractor = Bam_alignment_extractor(bam_file)
        pretty_alignments = bam_extractor.get_read_alignments(contig_acc, pretty=True)
        
        grouped_alignments = self._group_alignments_by_read_name(pretty_alignments)

        mp_counter = MultiPathCounter()
        
        for read_name in grouped_alignments:
            #print("{}\t{}".format(read_name, len(grouped_alignments[read_name])))
            paths_list = list()
            for pretty_alignment in grouped_alignments[read_name]:
                path = self._map_read_to_graph(pretty_alignment.get_pretty_alignment_segments())
                #print("pretty_alignment: {} maps to graph path: {}".format(pretty_alignment, path))
                if path and path != SPACER:
                    paths_list.append(path)

            if paths_list:
                mp = MultiPath(self._splice_graph, paths_list)
                #print("paths_list: {} -> mp: {}".format(paths_list, mp))
                mp_counter.add(mp)

        print(mp_counter)
            
        return
    
         

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
                if path_part:
                    path_part.extend(self._map_segment_to_graph_TERMINAL(segment))
            else:
                # internal segment
                path_part = self._get_intron_node_id(alignment_segments[i-1], segment)
                if path_part:
                    path_part.extend(self._map_segment_to_graph_INTERNAL(segment))

            #print("segment: {}  mapped to {}".format(segment, path_part))
                    
            if path_part:
                path.extend(path_part)
                #print("\tpath extended to: {}".format(path))
            else:
                if len(path) == 0 or path[-1] != SPACER:
                    path.append(SPACER) # spacer


        # trim any terminal spacer
        if path[-1] == SPACER:
            path = path[:-1]
                    
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
            # check for overlap and not extending beyond feature rend
            if (segment[0] < exon_segment._rend and
                segment[1] > exon_segment._lend and
                segment[0] >= exon_segment._lend and
                segment[1] <= exon_segment._rend):

                path.append(exon_segment.get_id())

        return path

    
