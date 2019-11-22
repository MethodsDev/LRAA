#!/usr/bin/env python
# encoding: utf-8

import sys, os, re
import subprocess
import logging
import string
import pysam

logger = logging.getLogger(__name__)


class splice_graph:

    # ---------------
    # class variables
    _read_aln_gap_merge_int = 10
    _max_contig_length = 1e10

    
    
    def __init__(self):

        # ------------------
        # instance variables
        self._genome_fasta_filename = ""
        self._alignments_bam_filename = ""
        
        self._contig_acc = ""
        self._contig_seq_str = ""
        self._contig_base_cov = []
    
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


        # get genome contig sequence
        contig_seq_str = self._retrieve_contig_seq()
        contig_len = len(contig_seq_str)
        logging.info("initing coverage array of len: {}".format(contig_len)) 
        if (contig_len > splice_graph._max_contig_length):
            raise RuntimeError("contig length {} exceeds maximum allowed {}".format(contig_len, splice_graph._max_contig_length))

        
        # init depth of coverage array
        self._contig_base_cov = [0 for i in range(0,contig_len)]
        
        samfile = pysam.AlignmentFile(alignments_bam_file, "rb")
        for read in samfile.fetch(contig_acc):
            alignment_segments = self._get_alignment_segments(read)

            if len(alignment_segments) > 1:
                # ensure proper consensus splice sites.
                if not self._matches_splicing_consensus(alignment_segments):
                    continue
                

            
            
        
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

        # print(contig_seq_str)
        
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
    
        
    def _matches_splicing_consensus(self, alignment_segments):

        genome_seq = self._contig_seq_str

        top_strand_agreement_count = 0
        bottom_strand_agreement_count = 0
        
        for i in range(len(alignment_segments) -1):
            seg_left_rend = alignment_segments[i][1]
            seg_right_lend = alignment_segments[i+1][0]
            
            dinuc_left = genome_seq[seg_left_rend] + genome_seq[seg_left_rend + 1]
                        
            dinuc_right = genome_seq[seg_right_lend - 2 ] + genome_seq[seg_right_lend - 1 ]

            dinuc_combo = dinuc_left + dinuc_right
            
            if dinuc_combo in self._splice_dinucs_top_strand:
                top_strand_agreement_count += 1
            elif dinuc_combo in self._splice_dinucs_bottom_strand:
                bottom_strand_agreement_count += 1
            else:
                return False


        if top_strand_agreement_count > 0 and bottom_strand_agreement_count > 0:
            # inconsistent orientation of splicing events
            return False
        elif top_strand_agreement_count > 0 or bottom_strand_agreement_count > 0:
            # all consistent splicing orientations
            return True
        else:
            raise RuntimeError("splicing analysis error... shouldn't happen")


    
        
