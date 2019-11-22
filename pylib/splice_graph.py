#!/usr/bin/env python
# encoding: utf-8

import sys, os, re
import subprocess
import logging
import string
import pysam
from collections import defaultdict

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

        pass
    
