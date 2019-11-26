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

logger = logging.getLogger(__name__)

DEBUG=True


class GenomeFeature:

    def __init__(self, lend, rend):
        self._lend = lend
        self._rend = rend
        self._id = "__id_not_set__"
        self._contig_acc = "__contig_acc_not_set__"
        return
        
    def get_coords(self):
        return(self._lend, self._rend)

    def get_read_support(self):
        # implement your own!
        return(0)

    def get_feature_length(self):
        return(self._rend - self._lend + 1)

    def get_bed_row(self, pad=0):
        return("\t".join([ str(x) for x in [self._contig_acc, self._lend - pad, self._rend + pad, self._id, self.get_read_support()] ])) 
    
        
class Intron(GenomeFeature):

    intron_id_counter = 0
    
    def __init__(self, lend, rend, orient, count):
        super().__init__(lend, rend)
        self._orient = orient
        self._count = count

        Intron.intron_id_counter += 1
        self._id = "I:{}".format(Intron.intron_id_counter)

        return

    def get_read_support(self):
        return self._count
    
    def __repr__(self):
        return("Intron: {} {}-{} count:{}".format(self._id, self._lend, self._rend, self._count))
        
    
class Exon(GenomeFeature):

    exon_id_counter = 0

    def __init__(self, lend, rend, mean_coverage):
        super().__init__(lend, rend)
        self._mean_coverage = mean_coverage

        Exon.exon_id_counter += 1
        self._id = "E:{}".format(Exon.exon_id_counter)

        return

    def get_read_support(self):
        return self._mean_coverage
    
    def __repr__(self):
        return("Exon: {} {}-{} mean_cov:{}".format(self._id, self._lend, self._rend, self._mean_coverage))

class splice_graph:

    # ---------------
    # class variables
    # ---------------
    
    _read_aln_gap_merge_int = 10
    _max_contig_length = 1e10

    # noise filtering params
    _min_alt_splice_freq = 0.05
    _max_intron_length_for_exon_segment_filtering = 10000
    
    
    def __init__(self):

        # ------------------
        # instance variables
        self._genome_fasta_filename = ""
        self._alignments_bam_filename = ""
        
        self._contig_acc = ""
        self._contig_seq_str = ""
        
        self._contig_base_cov = list()
        self._introns = defaultdict(int)
        
        self._splice_graph = None  # becomes networkx digraph

        self._splice_dinucs_top_strand = { "GTAG", "GCAG", "ATAC" }
        self._splice_dinucs_bottom_strand = {"CTAC", "CTGC", "GTAT" } # revcomp of top strand dinucs
        
        return


    def set_read_aln_gap_merge(self, read_aln_gap_merge_int):

        self._read_aln_gap_merge_int = read_aln_gap_merge_int

        return
    

    def build_splice_graph_for_contig(self, contig_acc, genome_fasta_file, alignments_bam_file):


        logger.info("creating splice graph for {} leveraging fasta {} and bam {}".format(contig_acc,
                                                                                         genome_fasta_file,
                                                                                         alignments_bam_file))
        self._contig_acc = contig_acc
        self._genome_fasta_filename = genome_fasta_file
        self._alignments_bam_filename = alignments_bam_file

        ## do the work:
        
        self._extract_contig_coverage()
        self._extract_introns()
        
        draft_splice_graph = self._build_draft_splice_graph()

        draft_splice_graph = self._prune_lowly_expressed_intron_overlapping_exon_segments(draft_splice_graph)

        draft_splice_graph = self._prune_disconnected_introns(draft_splice_graph)

        self._splice_graph = draft_splice_graph
        
        return self._splice_graph
    
    

    def _extract_contig_coverage(self):

        ## Contig Depth Array Capture
        # get genome contig sequence
        contig_seq_str = self._retrieve_contig_seq()
        contig_len = len(contig_seq_str)
        logging.info("initing coverage array of len: {}".format(contig_len)) 
        if (contig_len > splice_graph._max_contig_length):
            raise RuntimeError("contig length {} exceeds maximum allowed {}".format(contig_len, splice_graph._max_contig_length))

        
        # init depth of coverage array
        self._contig_base_cov = [0 for i in range(0,contig_len)]

        return

    
    def _extract_introns(self):
        
        ## Intron Capture
        # parse read alignments, capture introns and genome coverage info.
        samfile = pysam.AlignmentFile(self._alignments_bam_filename, "rb")
        for read in samfile.fetch(self._contig_acc):
            alignment_segments = self._get_alignment_segments(read)

            if len(alignment_segments) > 1:
                # ensure proper consensus splice sites.
                introns_list = self._get_introns_matching_splicing_consensus(alignment_segments)
                if introns_list is None:
                    continue

                for intron in introns_list:
                    self._introns[intron] += 1

                
            # add to coverage
            for segment in alignment_segments:
                for i in range(*segment):
                    self._contig_base_cov[i] += 1
        
        return

    

    def _retrieve_contig_seq(self):

        contig_seq_str = subprocess.check_output("samtools faidx {} {}".format(self._genome_fasta_filename,
                                                                               self._contig_acc),
                                                 shell=True,
                                                 encoding="utf-8")
        
        contig_seq_str = contig_seq_str.upper()
        contig_seq_str = contig_seq_str.split("\n")
        contig_seq_str = contig_seq_str[1:]
        contig_seq_str = "".join(contig_seq_str)
                
        contig_seq_str = re.sub("\s", "", contig_seq_str) # just in case
        
        self._contig_seq_str = contig_seq_str

        return(contig_seq_str)


    def _get_alignment_segments(self, pysam_read_alignment):
        
        aligned_pairs = pysam_read_alignment.get_blocks()

        ## merge adjacent blocks within range.
        alignment_segments = list()
        alignment_segments.append(list(aligned_pairs.pop(0)))

        for aligned_pair in aligned_pairs:
            # extend earlier stored segment or append new one
            if aligned_pair[0] - alignment_segments[-1][1] < splice_graph._read_aln_gap_merge_int:
                # extend rather than append
                alignment_segments[-1][1] = aligned_pair[1]
            else:
                # append, as too far apart from prev
                alignment_segments.append(list(aligned_pair))


        return alignment_segments
    
        
    def _get_introns_matching_splicing_consensus(self, alignment_segments):
        
        genome_seq = self._contig_seq_str

        top_strand_agreement_count = 0
        bottom_strand_agreement_count = 0

        introns_list = list()
        
        for i in range(len(alignment_segments) -1):
            seg_left_rend = alignment_segments[i][1]  # exon coord not inclusive
            seg_right_lend = alignment_segments[i+1][0] # exon coord inclusive

            intron_lend = seg_left_rend
            intron_rend = seg_right_lend - 1
            
            introns_list.append( (intron_lend, intron_rend) )
            
            dinuc_left = genome_seq[intron_lend] + genome_seq[intron_lend + 1]
                        
            dinuc_right = genome_seq[intron_rend - 1] + genome_seq[intron_rend]

            dinuc_combo = dinuc_left + dinuc_right
            #print(dinuc_combo)
            
            if dinuc_combo in self._splice_dinucs_top_strand:
                top_strand_agreement_count += 1
            elif dinuc_combo in self._splice_dinucs_bottom_strand:
                bottom_strand_agreement_count += 1
            else:
                return None


        if top_strand_agreement_count > 0 and bottom_strand_agreement_count > 0:
            # inconsistent orientation of splicing events
            return None
        elif top_strand_agreement_count > 0 or bottom_strand_agreement_count > 0:
            # all consistent splicing orientations
            return introns_list
        else:
            raise RuntimeError("splicing analysis error... shouldn't happen")


    
    def _build_draft_splice_graph(self):

        ## do some intron pruning
        self._prune_likely_false_introns()

        ## segment genome coverage
        exon_segments = self._segment_exon_by_coverage_n_splicing()

        draft_splice_graph = nx.DiGraph()

        ## add intron nodes.
        lend_to_intron = defaultdict(list)
        rend_to_intron = defaultdict(list)

        for intron in self._introns:
            intron_lend, intron_rend = intron
            count = self._introns[intron]

            intron_obj = Intron(intron_lend, intron_rend, "?", count)
            lend_to_intron[intron_lend].append(intron_obj)
            rend_to_intron[intron_rend].append(intron_obj)

            draft_splice_graph.add_node(intron_obj)


        ## add exon nodes and connect to introns.
        
        for exon in exon_segments:
            exon_lend, exon_rend = exon

            exon_mean_cov = self._get_mean_coverage(exon_lend, exon_rend)
            
            exon_obj = Exon(exon_lend, exon_rend, exon_mean_cov)

            draft_splice_graph.add_node(exon_obj)

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
        
        return draft_splice_graph            
        
        
        

    def _prune_likely_false_introns(self):

        self._prune_spurious_introns_shared_boundary("left")
        self._prune_spurious_introns_shared_boundary("right")


        
        
    def _prune_spurious_introns_shared_boundary(self, left_or_right):

        idx = 0 if left_or_right == "left" else 1
        
        introns_shared_coord = defaultdict(list)
        
        for intron in self._introns:
            intron_left = intron[idx]
            introns_shared_coord[intron_left].append(intron)

        introns_to_delete = set()
        
        # see if there are relatively poorly supported junctions
        for intron_list in introns_shared_coord.values():
            intron_list = sorted(intron_list, key=lambda x: self._introns[x])
            most_supported_intron = intron_list.pop()
            most_supported_intron_abundance = self._introns[most_supported_intron]
            for alt_intron in intron_list:
                alt_intron_abundance = self._introns[alt_intron]
                alt_intron_relative_freq = alt_intron_abundance / most_supported_intron_abundance
                #print("alt intron: {}".format(alt_intron) + " has rel freq: {}".format(alt_intron_relative_freq))
                if (alt_intron_relative_freq < splice_graph._min_alt_splice_freq):
                    introns_to_delete.add(alt_intron)

        logger.info("removing {} low frequency introns with shared {} coord".format(len(introns_to_delete), left_or_right))
        
        for intron in introns_to_delete:
            del self._introns[intron]
            
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
        if DEBUG:
            ofh = open("__splice_neighbor_cov_ratios.dat", 'w')
        
        for intron in self._introns:
            lend, rend = intron
            intron_abundance = self._introns[intron]
            
            A_mean_cov = self._get_mean_coverage(lend - coverage_window_length, lend - 1)
            B_mean_cov = self._get_mean_coverage(lend, lend + coverage_window_length)

            ratio_B_A = (B_mean_cov + pseudocount)  / (A_mean_cov + pseudocount)
                        
            C_mean_cov = self._get_mean_coverage(rend - coverage_window_length, rend)
            D_mean_cov = self._get_mean_coverage(rend + 1, rend + coverage_window_length)
            
            ratio_C_D = (C_mean_cov + pseudocount) / (D_mean_cov + pseudocount)

            if DEBUG:
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

        for intron in self._introns:
            lend, rend = intron
            left_splice_sites.add(lend)
            right_splice_sites.add(rend)


        exon_segments = list()
        exon_seg_start = None
        for i in range(0, len(self._contig_base_cov)):
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
                elif i in left_splice_sites:
                    exon_segments.append([exon_seg_start, i-1])
                    exon_seg_start = i
                elif i in right_splice_sites:
                    exon_segments.append([exon_seg_start, i])
                    exon_seg_start = None


        # get last one if it runs off to the end of the contig
        if exon_seg_start is not None:
            exon_segments.append([exon_seg_start, len(self._contig_base_cov)-1])


        if DEBUG:
            # write exon list to file
            with open("__exon_regions.init.bed", 'w') as ofh:
                for segment in exon_segments:
                    ofh.write("\t".join([self._contig_acc, str(segment[0]-1), str(segment[1]+1)]) + "\n")

            with open("__introns.init.bed", 'w') as ofh:
                for intron in self._introns:
                    intron_lend, intron_rend = intron
                    intron_support = self._introns[intron]
                    
                    ofh.write("\t".join([self._contig_acc, str(intron_lend), str(intron_rend), "{}-{}".format(intron_lend, intron_rend), str(intron_support)]) + "\n")
                
        return exon_segments
    

    def _prune_lowly_expressed_intron_overlapping_exon_segments(self, draft_splice_graph):

        intron_objs = list()
        exon_segment_objs = list()

        for node in draft_splice_graph:
            if type(node) == Intron:
                intron_objs.append(node)
            elif type(node) == Exon:
                exon_segment_objs.append(node)
            else:
                raise RuntimeError("Error, not identifying node: {} as Exon or Intron type - instead {} ".format(node, type(node)))

        #intron_objs = sorted(intron_objs, key=lambda x: x._lend)
        #exon_segment_objs = sorted(exon_segment_objs, lambda x: x._lend)

        # build interval tree for exon segments.

        exon_itree = itree.IntervalTree()
        for exon_seg in exon_segment_objs:
            exon_lend,exon_rend = exon_seg.get_coords()
            exon_itree[exon_lend:exon_rend+1] = exon_seg
        

        exons_to_purge = set()
            
        min_intron_cov_for_filtering = 1 / splice_graph._min_alt_splice_freq + 1
        
        for intron in intron_objs:

            if intron.get_read_support() < min_intron_cov_for_filtering:
                continue
            if intron.get_feature_length() > splice_graph._max_intron_length_for_exon_segment_filtering:
                continue
            
            
            intron_lend,intron_rend = intron.get_coords()
            overlapping_exon_segs = exon_itree[intron_lend:intron_rend+1]
            #print("Intron: {}".format(intron))
            
            for overlapping_exon_seg_iv in overlapping_exon_segs:
                #print("\toverlaps: {}".format(overlapping_exon_seg))

                overlapping_exon_seg = overlapping_exon_seg_iv.data
                
                if overlapping_exon_seg.get_read_support() / intron.get_read_support() < splice_graph._min_alt_splice_freq:
                    exons_to_purge.add(overlapping_exon_seg)


        logger.info("-removing {} lowly expressed exon segments based on intron overlap".format(len(exons_to_purge)))
        if exons_to_purge:
            draft_splice_graph.remove_nodes_from(exons_to_purge)
            

        return(draft_splice_graph)



    def _prune_disconnected_introns(self, draft_splice_graph):

        introns_to_remove = list()
        for node in draft_splice_graph:
            if type(node) == Intron:
                ## check it has at least one parent and one child
                if ( len(list(draft_splice_graph.predecessors(node))) == 0
                         or
                     len(list(draft_splice_graph.successors(node))) == 0 ):

                    introns_to_remove.append(node)

        logger.info("-pruning {} now disconnected introns".format(len(introns_to_remove)))

        draft_splice_graph.remove_nodes_from(introns_to_remove)

        return draft_splice_graph
    

    def write_intron_exon_splice_graph_bed_files(self, output_prefix):

        exons_bed_file = "{}.exons.bed".format(output_prefix)
        introns_bed_file = "{}.introns.bed".format(output_prefix)

        exons_ofh = open(exons_bed_file, 'w')
        introns_ofh = open(introns_bed_file, 'w')
        
        for node in self._splice_graph:
            if type(node) == Exon:
                exons_ofh.write(node.get_bed_row(pad=1) + "\n")
            elif type(node) == Intron:
                introns_ofh.write(node.get_bed_row(pad=1) + "\n")
            else:
                raise RuntimeError("not intron or exon object... bug... ")

        exons_ofh.close()
        introns_ofh.close()

        return
    
                
