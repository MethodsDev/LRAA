#!/usr/bin/env python3
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
import PASA_SALRAA_Globals

logger = logging.getLogger(__name__)


class Splice_graph:

    # ---------------
    # class variables
    # ---------------
    
    _read_aln_gap_merge_int = 10 ## internal alignment gap merging
    _inter_exon_segment_merge_dist = 50  ## unbranched exon segments within this range get merged into single segments.
    _max_genomic_contig_length = 1e10

    # noise filtering params
    _min_alt_splice_freq = 0.01
    _min_alt_unspliced_freq = 0.20
    _max_intron_length_for_exon_segment_filtering = 10000
    _min_intron_support = 1
    _min_terminal_splice_exon_anchor_length = 15


    _remove_unspliced_introns = False #True
    
    def __init__(self):

        # ------------------
        # instance variables
        #self._alignments_bam_filename = ""
        
        self._contig_acc = ""
        self._contig_seq_str = ""
        self._contig_len = 0
        
        self._contig_base_cov = list()
        
        self._splice_graph = None  # becomes networkx digraph TODO://rename this var as confusing when using _splice_graph as this obj in other modules
        
        self._node_id_to_node = dict()
        self._itree_exon_segments = itree.IntervalTree()
        self._intron_objs = dict() # "lend:rend" => intron_obj

        ## connected components (genes)
        self._components = list() # ordered list of graph components
        self._node_id_to_component = dict() # node_id -> component index
        
        
        return

    @classmethod
    def init_sg_params(cls, params):
        cls._read_aln_gap_merge_int = params['read_aln_gap_merge_int'] ## internal alignment gap merging
        cls._inter_exon_segment_merge_dist = params['inter_exon_segment_merge_dist']  ## unbranched exon segments within this range get merged into single segments.
        cls._max_genomic_contig_length = params['max_genomic_contig_length']

        # noise filtering params
        cls._min_alt_splice_freq = params['min_alt_splice_freq']
        cls._min_alt_unspliced_freq = params['min_alt_unspliced_freq']
        cls._max_intron_length_for_exon_segment_filtering = params['max_intron_length_for_exon_segment_filtering']
        cls._min_intron_support = params['min_intron_support']
        cls._min_terminal_splice_exon_anchor_length = params['min_terminal_splice_exon_anchor_length']

        cls._remove_unspliced_introns = params['remove_unspliced_introns']

    

    def get_contig_acc(self):
        return self._contig_acc

    
    def set_read_aln_gap_merge(self, read_aln_gap_merge_int):

        self._read_aln_gap_merge_int = read_aln_gap_merge_int

        return

    
    def get_intron_node_obj(self, intron_lend, intron_rend):

        intron_token = "{}:{}".format(intron_lend, intron_rend)
                
        if intron_token in self._intron_objs:
            return self._intron_objs[intron_token]

        return None


    def get_node_obj_via_id(self, node_id):

        if node_id not in self._node_id_to_node:
            raise RuntimeError("Error, node_id: {} not recognized in splice graph".format(node_id))
        
        return(self._node_id_to_node[node_id])


    def get_successors(self, node):
        return list(self._splice_graph.successors(node))

    def get_predecessors(self, node):
        return list(self._splice_graph.predecessors(node))
            
    
    def get_overlapping_exon_segments(self, range_lend, range_rend):

        overlapping_exon_segments = list()

        for overlapping_interval in self._itree_exon_segments[range_lend:range_rend]:
            overlapping_exon_segments.append(overlapping_interval.data)

        return overlapping_exon_segments
    

    def build_splice_graph_for_contig(self, contig_acc, contig_seq_str, alignments_bam_file, region_lend, region_rend, input_transcripts=None):

        logger.info("creating splice graph for {} leveraging bam {}".format(contig_acc,
                                                                            alignments_bam_file))
        
        self._contig_acc = contig_acc
        self._contig_seq_str = contig_seq_str
        self._contig_seq_len = len(contig_seq_str)
        #self._alignments_bam_filename = alignments_bam_file
        self._region_lend = region_lend
        self._region_rend = region_rend
        
        ## do the work:
        
        self._initialize_contig_coverage() # populates self._contig_base_cov

        ##---------------------------------------------------
        ## intron extraction and exonic base coverage defined.
        # -requires min intron-supporting reads
        # -excludes introns w/ heavily unbalanced splice site support (via Splice_graph._min_alt_splice_freq setting)
        # - stores intron objs in self._intron_objs
        # - base coverage incremented under self._contig_base_cov

        if alignments_bam_file is not None:
            self._populate_exon_coverage_and_extract_introns(alignments_bam_file)

        # incorporate guide structures if provided
        if input_transcripts:
            self._integrate_input_transcript_structures(input_transcripts)

        ##--------------------------------------------------------------------------------
        # initializes self._splice_graph
        # -defines exons by segmenting genomic coverage based on intron splice coordinates
        # - constructs self._splice_graph as nx.DiGraph()
        self._build_draft_splice_graph() 

        if self.is_empty():
            return None

        
        if PASA_SALRAA_Globals.DEBUG:
            self.write_intron_exon_splice_graph_bed_files("__prefilter", pad=0)
            self.describe_graph("__prefilter.graph")
                

        ##----------------------------------------------
        ## Refine the splice graph
        ## removes exon segments, not introns, also removes unspliced introns if Splice_graph._remove_unspliced_introns flag is set.

        self._prune_lowly_expressed_intron_overlapping_exon_segments()  
        
        self._merge_neighboring_proximal_unbranched_exon_segments()
        
        # self._prune_exon_spurs_at_introns() ## TODO:// implement this in a data-driven way.

        if Splice_graph._remove_unspliced_introns:
            self._prune_unspliced_introns()
        
        self._prune_disconnected_introns()
                
        self._finalize_splice_graph()

        if PASA_SALRAA_Globals.DEBUG:
            self.write_intron_exon_splice_graph_bed_files("__final_graph", pad=0)
            self.describe_graph("__final.graph")

        #print(self._itree_exon_segments )


        connected_components = list(nx.connected_components(self._splice_graph.to_undirected()))

        self._components = connected_components
        for i,component in enumerate(self._components):
            for node in component:
                id = node.get_id()
                logger.debug(f"assigning node {id} to component {i}")
                self._node_id_to_component[id] = i
        

        return self._splice_graph
    

    
    def _node_has_successors(self, node):

        if len(list(self._splice_graph.successors(node))) > 0:
            return True
        else:
            return False

    def _node_has_predecessors(self, node):

        if (len(list(self._splice_graph.predecessors(node)))) > 0:
            return True
        else:
            return False
        

    def _get_exon_and_intron_nodes(self):

        intron_objs = list()
        exon_segment_objs = list()

        for node in self._splice_graph:
            if type(node) == Intron:
                intron_objs.append(node)
            elif type(node) == Exon:
                exon_segment_objs.append(node)
            else:
                raise RuntimeError("Error, not identifying node: {} as Exon or Intron type - instead {} ".format(node, type(node)))

        exon_segment_objs = sorted(exon_segment_objs, key=lambda x: x._lend)
        intron_objs = sorted(intron_objs, key=lambda x: x._lend)

        return (exon_segment_objs, intron_objs)

    

    def _initialize_contig_coverage(self):

        ## Contig Depth Array Capture
        # get genome contig sequence
        contig_seq_str = self._contig_seq_str
        contig_len = self._contig_seq_len
        logging.info("initing coverage array of len: {}".format(contig_len)) 
        if (contig_len > Splice_graph._max_genomic_contig_length):
            raise RuntimeError("genomic contig length {} exceeds maximum allowed {}".format(contig_len, Splice_graph._max_genomic_contig_length))

        
        # init depth of coverage array
        self._contig_base_cov = [0 for i in range(0,contig_len+1)]

        return

    
    def _populate_exon_coverage_and_extract_introns(self, bam_filename):
        
        ## Intron Capture

        intron_counter = defaultdict(int)

        intron_to_read_types = defaultdict(set)
        
        intron_splice_site_support = defaultdict(int)
        
        bam_extractor = Bam_alignment_extractor(bam_filename)
        
        # get read alignments
        # - illumina and pacbio reads filtered based on tech-specific min per_id
        # - pretty alignments: store the pysam alignment record along with inferred transcript exons segments.
        pretty_alignments = bam_extractor.get_read_alignments(self._contig_acc,
                                                              region_lend=self._region_lend, region_rend=self._region_rend,
                                                              pretty=True)
        logger.info("-got {} pretty alignments.".format(len(pretty_alignments)))
        
        total_read_alignments_used = 0
        
        for pretty_alignment in pretty_alignments:

            alignment_segments = pretty_alignment.get_pretty_alignment_segments()
            #print("Pretty alignment segments: " + str(alignment_segments))

            read_type = pretty_alignment.get_read_type()
            
            if len(alignment_segments) > 1:
                # ensure proper consensus splice sites.
                introns_list = self._get_introns_matching_splicing_consensus(alignment_segments)
                #print("introns list: " + str(introns_list))
                for intron in introns_list:
                    intron_counter[intron] += 1
                    intron_to_read_types[intron].add(read_type)
                    
                    intron_lend,intron_rend,splice_orient = intron
                    intron_splice_site_support[intron_lend] += 1
                    intron_splice_site_support[intron_rend] += 1

            total_read_alignments_used += 1
            # add to coverage
            for segment in alignment_segments:
                for i in range(segment[0], segment[1] + 1):
                    if i > self._contig_seq_len:
                        break
                    
                    self._contig_base_cov[i] += 1

        logger.info("-total read alignments used: {}".format(total_read_alignments_used))
        
        # retain only those introns that meet the min threshold
        for intron_coords, count in intron_counter.items():

            read_types = intron_to_read_types[intron_coords]
            
            if count >= Splice_graph._min_intron_support:
                ## check splice support
                intron_lend, intron_rend, intron_orient = intron_coords
                splice_support_left = intron_splice_site_support[intron_lend]
                splice_support_right = intron_splice_site_support[intron_rend]

                min_support = min(splice_support_left, splice_support_right)
                max_support = max(splice_support_left, splice_support_right)

                if min_support/max_support >= Splice_graph._min_alt_splice_freq:
                    intron_obj = Intron(self._contig_acc, intron_lend, intron_rend, intron_orient, count)
                    intron_obj.add_read_types(list(read_types))
                    
                    intron_coords_key = "{}:{}".format(intron_lend, intron_rend)
                    self._intron_objs[intron_coords_key] = intron_obj
            
        
                
        return


    
    def _integrate_input_transcript_structures(self, transcripts):

        """
        Fold in the reference annotations:
        - add introns where they're missing
        - add base coverage where it's missing.
        
        """
        
                
        ## Exon Coverage and Intron Capture

        for transcript in transcripts:
            orient = transcript.get_orient()
            exon_segments = transcript.get_exon_segments()
            last_rend = None

            # add coverage for exonic region
            for exon_segment in exon_segments:
                lend, rend = exon_segment
                for i in range(lend, rend + 1):
                    if i > self._contig_seq_len:
                        break
                    if self._contig_base_cov[i] < 1:
                        self._contig_base_cov[i] = 1

                # add missing introns.:
                if last_rend is not None:
                    intron_lend = last_rend + 1
                    intron_rend = lend - 1
                    intron_coords_key = "{}:{}".format(intron_lend, intron_rend)

                    if intron_coords_key not in self._intron_objs:
                        intron_obj = Intron(self._contig_acc, intron_lend, intron_rend, orient, 1)
                        intron_obj.add_read_types(['input_transcript']) ## TODO:// is this necessary? if not, remove.
                        self._intron_objs[intron_coords_key] = intron_obj


                last_rend = rend

        
        return



    
    def _get_introns_matching_splicing_consensus(self, alignment_segments):
        
        genome_seq = self._contig_seq_str

        top_strand_agreement_count = 0
        bottom_strand_agreement_count = 0

        introns_list = list()
        
        for i in range(len(alignment_segments) -1):
            seg_left_rend = alignment_segments[i][1]  # exon coord not inclusive
            seg_right_lend = alignment_segments[i+1][0] # exon coord inclusive

            intron_lend = seg_left_rend + 1
            intron_rend = seg_right_lend - 1

            splice_type = Intron.check_canonical_splicing(intron_lend, intron_rend, genome_seq)

            if splice_type is not None:
                introns_list.append( (intron_lend, intron_rend, splice_type) )

        #
        #    if splice_type is None:
        #        continue
        #    elif splice_type == '+':
        #        top_strand_agreement_count += 1
        #    elif splice_type == '-':
        #        bottom_strand_agreement_count += 1
        #    else:
        #        raise RuntimeError("not sure what splice type we have here...")
        #    
        # 
        #if top_strand_agreement_count > 0 and bottom_strand_agreement_count > 0:
        #    # inconsistent orientation of splicing events
        #    return None
        #elif top_strand_agreement_count > 0 or bottom_strand_agreement_count > 0:
        #    # all consistent splicing orientations
        #    return introns_list
        #else:
        #    raise RuntimeError("splicing analysis error... shouldn't happen")

        return introns_list

    
    def _build_draft_splice_graph(self):

        ## do some intron pruning
        self._prune_likely_false_introns()

        ## segment genome coverage
        exon_segments = self._segment_exon_by_coverage_n_splicing()

        draft_splice_graph = nx.DiGraph()

        ## add intron nodes.
        lend_to_intron = defaultdict(list)
        rend_to_intron = defaultdict(list)

        for intron in self._intron_objs.values():
            intron_lend, intron_rend = intron.get_coords()
            count = intron.get_read_support()

            lend_to_intron[intron_lend].append(intron)
            rend_to_intron[intron_rend].append(intron)

            draft_splice_graph.add_node(intron)


        ## add exon nodes and connect to introns.
        ## also connect adjacent exon segments.

        prev_exon_obj = None
        for exon in exon_segments:
            exon_lend, exon_rend = exon

            exon_mean_cov = self._get_mean_coverage(exon_lend, exon_rend)
            
            exon_obj = Exon(self._contig_acc, exon_lend, exon_rend, exon_mean_cov)

            draft_splice_graph.add_node(exon_obj)
            
            # connect adjacent exon segments.
            if prev_exon_obj is not None:
                if prev_exon_obj._rend + 1 == exon_obj._lend:
                    # adjacent segments.
                    draft_splice_graph.add_edge(prev_exon_obj, exon_obj)

            prev_exon_obj = exon_obj

            
            ## connect to introns where appropriate
            candidate_splice_left = exon_lend - 1
            if candidate_splice_left in rend_to_intron:
                introns = rend_to_intron[candidate_splice_left]
                for intron_obj in introns:
                    draft_splice_graph.add_edge(intron_obj, exon_obj)

            candidate_splice_right = exon_rend + 1
            if candidate_splice_right in lend_to_intron:
                introns = lend_to_intron[candidate_splice_right]
                for intron_obj in lend_to_intron[candidate_splice_right]:
                    draft_splice_graph.add_edge(exon_obj, intron_obj)


        self._splice_graph = draft_splice_graph

        return
        
        
        

    def _prune_likely_false_introns(self):

        self._prune_spurious_introns_shared_boundary("left")
        self._prune_spurious_introns_shared_boundary("right")


        
        
    def _prune_spurious_introns_shared_boundary(self, left_or_right):

        idx = 0 if left_or_right == "left" else 1
        
        introns_shared_coord = defaultdict(list)
        
        for intron in self._intron_objs.values():
            intron_left = intron.get_coords()[idx]
            introns_shared_coord[intron_left].append(intron)

        introns_to_delete = set()
        
        # see if there are relatively poorly supported junctions
        for intron_list in introns_shared_coord.values():
            intron_list = sorted(intron_list, key=lambda x: x.get_read_support())
            most_supported_intron = intron_list.pop()
            most_supported_intron_abundance = most_supported_intron.get_read_support()
            for alt_intron in intron_list:
                alt_intron_abundance = alt_intron.get_read_support()
                alt_intron_relative_freq = alt_intron_abundance / most_supported_intron_abundance

                if (alt_intron_relative_freq < Splice_graph._min_alt_splice_freq):
                    logger.debug("alt intron: {}".format(alt_intron) + " has rel freq: {}".format(alt_intron_relative_freq))
                    introns_to_delete.add(alt_intron)

        logger.info("removing {} low frequency introns with shared {} coord".format(len(introns_to_delete), left_or_right))
        
        for intron in introns_to_delete:
            intron_coords = intron.get_coords()
            intron_key = "{}:{}".format(intron_coords[0], intron_coords[1])
            intron_obj = self._intron_objs[intron_key]
            if intron_obj.has_read_type("PBLR"):
                logger.debug("-retaining intron {} as having long read support".format(str(intron_obj)))
            else:
                logger.debug("removing intron: {} {}".format(intron_key, self._intron_objs[intron_key]))
                del self._intron_objs[intron_key]
            
        return

    
    def _prune_weak_splice_neighboring_segments(self):

        #                    I
        #              _____________
        #             /             \
        #   *********/****       ****\************
        #       A       B         C      D
        #

        coverage_window_length = 10  # maybe make this configurable.
        pseudocount = 1

        
        introns_to_delete = set()

        ofh = None
        if PASA_SALRAA_Globals.DEBUG:
            ofh = open("__splice_neighbor_cov_ratios.dat", 'w')
        
        for intron in self._intron_objs.values():
            lend, rend = intron.get_coords()
            intron_abundance = intron.get_read_support()
            
            A_mean_cov = self._get_mean_coverage(lend - coverage_window_length, lend - 1)
            B_mean_cov = self._get_mean_coverage(lend, lend + coverage_window_length)

            ratio_B_A = (B_mean_cov + pseudocount)  / (A_mean_cov + pseudocount)
                        
            C_mean_cov = self._get_mean_coverage(rend - coverage_window_length, rend)
            D_mean_cov = self._get_mean_coverage(rend + 1, rend + coverage_window_length)
            
            ratio_C_D = (C_mean_cov + pseudocount) / (D_mean_cov + pseudocount)

            if PASA_SALRAA_Globals.DEBUG:
                ofh.write("{}".format(lend) +
                          "\tI:{}".format(intron_abundance) +
                          "\tA:{:.3f}".format(A_mean_cov) +
                          "\tB:{:.3f}".format(B_mean_cov) +
                          "\tB/A:{:.3f}".format(ratio_B_A) +
                          
                          "\n{}".format(rend) +
                          "\tC:{:.3f}".format(C_mean_cov) +
                          "\tD:{:.3f}".format(D_mean_cov) +
                          "\tC/D:{:.3f}".format(ratio_C_D) )
            
        return


    
    def _get_mean_coverage(self, start, end):

        cov_sum = 0
        cov_len = 0
        
        for i in range(start, end+1):
            if i>= 0 and i < len(self._contig_base_cov):
                cov_len += 1
                cov_sum += self._contig_base_cov[i]

        if cov_len > 0:
            mean_cov = cov_sum / cov_len
            return mean_cov
        else:
            logger.warning("coverage requested to position {} extends beyond length of contig {}".format(end, self._contig_acc))
            return 0
        


    def _segment_exon_by_coverage_n_splicing(self):


        left_splice_sites = set()
        right_splice_sites = set()

        for intron in self._intron_objs.values():
            lend, rend = intron.get_coords()
            left_splice_sites.add(lend)
            right_splice_sites.add(rend)


        exon_segments = list()
        exon_seg_start = None
        for i in range(1, len(self._contig_base_cov)):
            if exon_seg_start is None:
                if self._contig_base_cov[i]:
                    # start new exon segment
                    exon_seg_start = i
            else:
                # handle current segment
                if self._contig_base_cov[i] == 0:
                    # stop segment, add to seg list
                    exon_segments.append([exon_seg_start, i-1])
                    exon_seg_start = None

            # splice breakpoint logic
            if i+1 in left_splice_sites:
                if exon_seg_start is not None:
                    exon_segments.append([exon_seg_start, i])
                exon_seg_start = None

            if i in right_splice_sites:
                if exon_seg_start is not None:
                    exon_segments.append([exon_seg_start, i])
                exon_seg_start = None

                
        # get last one if it runs off to the end of the contig
        if exon_seg_start is not None:
            exon_segments.append([exon_seg_start, len(self._contig_base_cov)])


        if PASA_SALRAA_Globals.DEBUG:
            # write exon list to file
            with open("__exon_regions.init.bed", 'w') as ofh:
                for segment in exon_segments:
                    ofh.write("\t".join([self._contig_acc, str(segment[0]-1), str(segment[1]+1)]) + "\n")

            with open("__introns.init.bed", 'w') as ofh:
                for intron in self._intron_objs.values():
                    intron_lend, intron_rend = intron.get_coords()
                    intron_support = intron.get_read_support()
                    
                    ofh.write("\t".join([self._contig_acc, str(intron_lend), str(intron_rend), "{}-{}".format(intron_lend, intron_rend), str(intron_support)]) + "\n")
                
        return exon_segments
    

    def _prune_lowly_expressed_intron_overlapping_exon_segments(self):

        draft_splice_graph = self._splice_graph
        
        intron_objs = list()
        exon_segment_objs = list()

        for node in draft_splice_graph:
            if type(node) == Intron:
                intron_objs.append(node)
            elif type(node) == Exon:
                exon_segment_objs.append(node)
            else:
                raise RuntimeError("Error, not identifying node: {} as Exon or Intron type - instead {} ".format(node, type(node)))

        # store splice junction coordinates.
        splice_coords = set()
        for intron_obj in intron_objs:
            lend, rend = intron_obj.get_coords()
            splice_coords.add(lend)
            splice_coords.add(rend)
        
            
        # build interval tree for exon segments.

        exon_itree = itree.IntervalTree()
        for exon_seg in exon_segment_objs:
            exon_lend,exon_rend = exon_seg.get_coords()
            exon_itree[exon_lend:exon_rend+1] = exon_seg
        

        exons_to_purge = set()
            
        min_intron_cov_for_filtering = 1 / Splice_graph._min_alt_unspliced_freq + 1
        
        for intron in intron_objs:

            if intron.get_read_support() < min_intron_cov_for_filtering:
                continue
            if intron.get_feature_length() > Splice_graph._max_intron_length_for_exon_segment_filtering:
                continue
            
            
            intron_lend,intron_rend = intron.get_coords()
            overlapping_exon_segs = exon_itree[intron_lend:intron_rend+1]
            #print("Intron: {}".format(intron))
            
            for overlapping_exon_seg_iv in overlapping_exon_segs:
                #print("\toverlaps: {}".format(overlapping_exon_seg))

                overlapping_exon_seg = overlapping_exon_seg_iv.data

                # see if connected by an intron
                exon_lend, exon_rend = overlapping_exon_seg.get_coords()
                if exon_lend - 1 in splice_coords or exon_rend + 1 in splice_coords:
                    continue
                

                if float(overlapping_exon_seg.get_read_support()) / float(intron.get_read_support()) < Splice_graph._min_alt_unspliced_freq:
                    logger.debug("-pruning {} as has exon_read_support:{} / intron_read_support:{} < min_alt_unspliced_freq: {}".format(overlapping_exon_seg,
                                                                                                                                        overlapping_exon_seg.get_read_support(),
                                                                                                                                        intron.get_read_support(),
                                                                                                                                        Splice_graph._min_alt_unspliced_freq) )
                    exons_to_purge.add(overlapping_exon_seg)
                    

        
        logger.info("-removing {} lowly expressed exon segments based on intron overlap".format(len(exons_to_purge)))
        if exons_to_purge:
            draft_splice_graph.remove_nodes_from(exons_to_purge)

        if PASA_SALRAA_Globals.DEBUG:
            exons_to_purge = list(exons_to_purge)
            exons_to_purge = sorted(exons_to_purge, key=lambda x: x._lend)
            with open("__exon_segments_to_purge.bed", 'w') as ofh:
                for exon in exons_to_purge:
                    ofh.write(exon.get_bed_row(pad=1) + "\n")

        return



    def _prune_disconnected_introns(self):

        draft_splice_graph = self._splice_graph
        
        introns_to_remove = list()
        for node in draft_splice_graph:
            if type(node) == Intron:
                ## check it has at least one parent and one child
                if ( len(list(draft_splice_graph.predecessors(node))) == 0
                         or
                     len(list(draft_splice_graph.successors(node))) == 0 ):

                    introns_to_remove.append(node)

        logger.info("-pruning {} now disconnected introns".format(len(introns_to_remove)))

        if PASA_SALRAA_Globals.DEBUG:
            with open("__pruned_disconnected_introns.bed", 'w') as ofh:
                for intron in introns_to_remove:
                    ofh.write(intron.get_bed_row(pad=1) + "\n")

        # remove introns
        for intron_to_remove in introns_to_remove:
            intron_lend,intron_rend = intron_to_remove.get_coords()
            intron_token = "{}:{}".format(intron_lend,intron_rend)
            del self._intron_objs[intron_token]
        
        draft_splice_graph.remove_nodes_from(introns_to_remove)

        return
    

    def write_intron_exon_splice_graph_bed_files(self, output_prefix, pad=0):

        exons_bed_file = "{}.exons.bed".format(output_prefix)
        introns_bed_file = "{}.introns.bed".format(output_prefix)

        exons_ofh = open(exons_bed_file, 'w')
        introns_ofh = open(introns_bed_file, 'w')

        exons_list = list()
        introns_list = list()
        
        for node in self._splice_graph:
            if type(node) == Exon:
                exons_list.append(node)
            elif type(node) == Intron:
                introns_list.append(node)
            else:
                raise RuntimeError("not intron or exon object... bug... ")


        exons_list = sorted(exons_list, key=lambda x: x._lend)
        introns_list = sorted(introns_list, key=lambda x: x._lend)

        for exon in exons_list:
            exons_ofh.write(exon.get_bed_row(pad=pad) + "\n")

        for intron in introns_list:
            introns_ofh.write(intron.get_bed_row(pad=pad) + "\n")
        
        exons_ofh.close()
        introns_ofh.close()

        return
    

    def describe_graph(self, outputfilename):

        ofh = open(outputfilename, 'w')
        
        nodes = list(self._splice_graph.nodes)

        nodes = sorted(nodes, key=lambda x: x._lend)

        for node in nodes:
            node_descr = ""
            preds = list(self._splice_graph.predecessors(node))
            if preds:
                pred_strs = list()
                for pred in preds:
                    #print(pred)
                    pred_strs.append(str(pred))
                node_descr += ">;<".join(pred_strs)
            else:
                node_descr += "."

            node_descr += "\t<" + str(node) + ">\t"

            succs = list(self._splice_graph.successors(node))
            if succs:
                succs_strs = list()
                for succ in succs:
                    succs_strs.append(str(succ))
                node_descr += ">;<".join(succs_strs)
            else:
                node_descr += "."


            ofh.write(node_descr + "\n")


        ofh.close()


        return

    
    def _merge_neighboring_proximal_unbranched_exon_segments(self):


        merged_node_ids = list()
        
        exon_segment_objs, intron_objs = self._get_exon_and_intron_nodes()
        
        prev_node = exon_segment_objs[0]

        for i in range(1, len(exon_segment_objs)):
            next_node = exon_segment_objs[i]

            if ( (not self._node_has_successors(prev_node))
                and
                (not self._node_has_predecessors(next_node))
                and next_node._lend - prev_node._rend - 1 < Splice_graph._inter_exon_segment_merge_dist):

                # merge next node into the prev node
                prev_node._rend = next_node._rend
                prev_node._mean_coverage = self._get_mean_coverage(prev_node._lend, prev_node._rend)

                for next_node_successor in self._splice_graph.successors(next_node):
                    self._splice_graph.add_edge(prev_node, next_node_successor)

                # remove next node
                merged_node_ids.append(next_node.get_id())
                self._splice_graph.remove_node(next_node)
            else:
                prev_node = next_node
                
        if PASA_SALRAA_Globals.DEBUG:
            with open("__merged_exons.list", "wt") as ofh:
                print("\n".join(merged_node_ids), file=ofh)

        
                
    def _prune_exon_spurs_at_introns(self):

        exon_segment_objs, intron_objs = self._get_exon_and_intron_nodes()


        #######           ---------------------
        #                /       ^intron^      \
        #    -----------/===                 ===\-------------------
        #                R_spur              L_spur


        ## //TODO: incorporate read coverage checks 
        
        def is_R_spur(exon_node):

            has_intron_successor = False
            for successor in self._splice_graph.successors(exon_node):
                if type(successor) == Intron:
                    has_intron_successor = True


            has_intron_predecessor = False
            has_alt_intron = False
            
            
            for predecessor in self._splice_graph.predecessors(exon_node):
                if type(predecessor) == Intron:
                    has_intron_predecessor = True
                
                elif type(predecessor) == Exon:
                    for successor in self._splice_graph.successors(predecessor):
                        if type(successor) == Intron:
                            has_alt_intron = True

            return (not has_intron_successor) and (not has_intron_predecessor) and has_alt_intron

            
        def is_L_spur(exon_node):

            has_intron_predecessor = False

            for predecessor in self._splice_graph.predecessors(exon_node):
                if type(predecessor) == Intron:
                    has_intron_predecessor = True
            

            has_intron_successor = False
            has_alt_intron = False

            for successor in self._splice_graph.successors(exon_node):
                if type(successor) == Intron:
                    has_intron_successor = True
                
                elif type(successor) == Exon:
                    for predecessor in self._splice_graph.predecessors(successor):
                        if type(predecessor) == Intron:
                            has_alt_intron = True

            return (not has_intron_predecessor) and (not has_intron_successor) and has_alt_intron
        

        exons_to_prune = list()
                
        for exon in exon_segment_objs:
            if exon.get_feature_length() >= Splice_graph._min_terminal_splice_exon_anchor_length:
                # long enough, we'll keep it for now.
                continue

            # ok shortie, must see if it's a spur
            if is_L_spur(exon) or is_R_spur(exon):
                exons_to_prune.append(exon)
                
        
        if exons_to_prune:
            logger.info("-removing {} exon spurs".format(len(exons_to_prune)))

            if PASA_SALRAA_Globals.DEBUG:
                with open("__pruned_exon_spurs.list", "wt") as ofh:
                    for exon in exons_to_prune:
                        print(exon.get_id(), file=ofh)

            self._splice_graph.remove_nodes_from(exons_to_prune)

        return


    def _finalize_splice_graph(self):
        
        ## store node ID to node object
        for node in self._splice_graph:
            self._node_id_to_node[ node.get_id() ] = node
            if type(node) == Exon:
                # store exon segments in itree for overlap queries
                lend, rend = node.get_coords()
                self._itree_exon_segments[lend:rend+1] = node
            
        return


    def _is_unspliced_exon_segment_artifact(self, exon, intron):


        if exon.get_read_support() < Splice_graph._min_alt_unspliced_freq * intron.get_read_support():
            return True
        
            
        intron_lend, intron_rend = intron.get_coords()
        exon_lend, exon_rend = exon.get_coords()

        
        # prune singleton exon segments falling in sufficiently expressed introns:
        if (not self._node_has_predecessors(exon)) and (not self._node_has_successors(exon)):
            return True

                

        return False



    def _prune_unspliced_introns(self):

        draft_splice_graph = self._splice_graph
        
        intron_objs = list()
        exon_segment_objs = list()

        for node in draft_splice_graph:
            if type(node) == Intron:
                intron_objs.append(node)
            elif type(node) == Exon:
                exon_segment_objs.append(node)
            else:
                raise RuntimeError("Error, not identifying node: {} as Exon or Intron type - instead {} ".format(node, type(node)))

        # build interval tree for exon segments.

        exon_itree = itree.IntervalTree()
        for exon_seg in exon_segment_objs:
            exon_lend,exon_rend = exon_seg.get_coords()
            exon_itree[exon_lend:exon_rend+1] = exon_seg
        

        exons_to_purge = set()

        ## should we restrict to certain introns here? min coverage? YES!!!
        
        for intron in intron_objs:
            intron_lend,intron_rend = intron.get_coords()
            overlapping_exon_segs = exon_itree[intron_lend:intron_rend+1]

            for overlapping_exon_seg in overlapping_exon_segs:
                exon_seg = overlapping_exon_seg.data
                if self._is_unspliced_exon_segment_artifact(exon=exon_seg, intron=intron):
                    exons_to_purge.add(exon_seg)

        logger.info("-removing {} likely unspliced exon segments based on intron overlap".format(len(exons_to_purge)))
        if exons_to_purge:
            draft_splice_graph.remove_nodes_from(exons_to_purge)

        if PASA_SALRAA_Globals.DEBUG:
            exons_to_purge = list(exons_to_purge)
            exons_to_purge = sorted(exons_to_purge, key=lambda x: x._lend)
            with open("__exon_segments_to_purge.bed", 'a') as ofh: # file should be already created based on earlier low expressed exon segments overlapping introns removal step
                for exon in exons_to_purge:
                    ofh.write(exon.get_bed_row(pad=1) + "\n")

        return


    def is_empty(self):
        return len(self._splice_graph) == 0

