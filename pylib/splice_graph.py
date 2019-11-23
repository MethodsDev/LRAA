#!/usr/bin/env python
# encoding: utf-8

import sys, os, re
import subprocess
import logging
import string
import pysam
from collections import defaultdict
import networkx as nx


logger = logging.getLogger(__name__)

DEBUG=True


class GenomeFeature:

    def __init__(self, lend, rend):
        self._lend = lend
        self._rend = rend
        return
        
    def get_coords(self):
        return(self._lend, self._rend)

    
class Intron(GenomeFeature):

    intron_id_counter = 0
    
    def __init__(self, lend, rend, orient, count):
        super().__init__(lend, rend)
        self._orient = orient
        self._count = count

        Intron.intron_id_counter += 1
        self._id = "I:{}".format(Intron.intron_id_counter)

        return

class Exon(GenomeFeature):

    exon_id_counter = 0

    def __init__(self, lend, rend, mean_coverage):
        super().__init__(lend, rend)
        self._mean_coverage = mean_coverage

        Exon.exon_id_counter += 1
        self._id = "E:{}".format(Exon.exon_id_counter)

        return



class splice_graph:

    # ---------------
    # class variables
    _read_aln_gap_merge_int = 10
    _max_contig_length = 1e10
    _min_alt_splice_freq = 0.05
    
    
    def __init__(self):

        # ------------------
        # instance variables
        self._genome_fasta_filename = ""
        self._alignments_bam_filename = ""
        
        self._contig_acc = ""
        self._contig_seq_str = ""
        
        self._contig_base_cov = list()
        self._introns = defaultdict(int)
        
        self._splice_graph = dict()

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
        
        self._extract_introns_and_contig_coverage()
        
        self._build_draft_splice_graph()

        
        return

    

    def _extract_introns_and_contig_coverage(self):
        
        # get genome contig sequence
        contig_seq_str = self._retrieve_contig_seq()
        contig_len = len(contig_seq_str)
        logging.info("initing coverage array of len: {}".format(contig_len)) 
        if (contig_len > splice_graph._max_contig_length):
            raise RuntimeError("contig length {} exceeds maximum allowed {}".format(contig_len, splice_graph._max_contig_length))

        
        # init depth of coverage array
        self._contig_base_cov = [0 for i in range(0,contig_len)]

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
        exome_segments = self._segment_exome_by_coverage_n_splicing()

        draft_splice_graph = nx.Graph()

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
        
        for exon in exome_segments:
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
        


    def _segment_exome_by_coverage_n_splicing(self):


        left_splice_sites = set()
        right_splice_sites = set()

        for intron in self._introns:
            lend, rend = intron
            left_splice_sites.add(lend)
            right_splice_sites.add(rend)


        exome_segments = list()
        exon_seg_start = None
        for i in range(0, len(self._contig_base_cov)):
            if exon_seg_start is None:
                if self._contig_base_cov[i]:
                    # start new exome segment
                    exon_seg_start = i
            else:
                # handle current segment
                if self._contig_base_cov[i] == 0:
                    # stop segment, add to seg list
                    exome_segments.append([exon_seg_start, i-1])
                    exon_seg_start = None
                elif i in left_splice_sites:
                    exome_segments.append([exon_seg_start, i-1])
                    exon_seg_start = i
                elif i in right_splice_sites:
                    exome_segments.append([exon_seg_start, i])
                    exon_seg_start = None


        # get last one if it runs off to the end of the contig
        if exon_seg_start is not None:
            exome_segments.append([exon_seg_start, len(self._contig_base_cov)-1])


        if DEBUG:
            # write exon list to file
            with open("__exome_regions.bed", 'w') as ofh:
                for segment in exome_segments:
                    ofh.write("\t".join([self._contig_acc, str(segment[0]-1), str(segment[1]+1)]) + "\n")

        return exome_segments
    
