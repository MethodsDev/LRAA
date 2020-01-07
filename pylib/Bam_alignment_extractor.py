#!/usr/bin/env python
# encoding: utf-8

import sys, os, re
import subprocess
import logging
import string
import pysam
from collections import defaultdict
from Pretty_alignment import Pretty_alignment

logger = logging.getLogger(__name__)

class Bam_alignment_extractor:

    # ---------------
    # class variables
    # ---------------
    
        
    def __init__(self, alignments_bam_filename):

        # ------------------
        # instance variables

        self._read_aln_gap_merge_int = 10 ## internal alignment gap merging
        self._min_terminal_splice_exon_anchor_length = 15
        self._min_read_aln_per_id = 98

        self._alignments_bam_filename = alignments_bam_filename

        self._pysam_reader = pysam.AlignmentFile(self._alignments_bam_filename, "rb")

        return


    def set_read_aln_gap_merge(self, read_aln_gap_merge_int):

        self._read_aln_gap_merge_int = read_aln_gap_merge_int

        return
    
    
    def get_read_alignments(self, contig_acc, pretty=False):        

        discarded_read_counter = defaultdict(int)

        read_alignments = list()
        
        # parse read alignments, capture introns and genome coverage info.

        for read in self._pysam_reader.fetch(contig_acc):
            
            if read.is_paired and not read.is_proper_pair:
                discarded_read_counter["improper_pair"] += 1
                continue

            if read.is_duplicate:
                discarded_read_counter["duplicate"] += 1
                continue

            if read.is_qcfail:
                discarded_read_counter["qcfail"] += 1
                continue

            if read.is_secondary:
                discarded_read_counter["secondary"] += 1
                continue
            
            # check read alignment percent identity
            cigar_stats = read.get_cigar_stats()
            aligned_base_count = cigar_stats[0][0]
            mismatch_count = None
            if read.has_tag("NM"):
                mismatch_count = int(read.get_tag("NM"))
            elif read.has_tag("nM"):
                mismatch_count = int(read.get_tag("nM"))
            if mismatch_count is not None:
                per_id = 100 - (mismatch_count/aligned_base_count)*100
                if per_id < self._min_read_aln_per_id:
                    discarded_read_counter["low_perID"] += 1
                    continue

            read_alignments.append(read)


        logger.info("reads discarded: {}".format(discarded_read_counter))
        
        if pretty:
            return self.get_pretty_alignments(read_alignments)
        else:
            return read_alignments


    def get_pretty_alignments(self, read_alignments_list):

        pretty_alignments = list()
        
        for read_alignment in read_alignments_list:
                    
            alignment_segments = self._get_alignment_segments(read_alignment)
            this_pretty_alignment = Pretty_alignment(read_alignment, alignment_segments)

            pretty_alignments.append(this_pretty_alignment)

        return pretty_alignments
              

    def _get_alignment_segments(self, pysam_read_alignment):
        
        aligned_pairs = pysam_read_alignment.get_blocks()

        #print(aligned_pairs)
        
        ## merge adjacent blocks within range.
        alignment_segments = list()
        alignment_segments.append(list(aligned_pairs.pop(0)))
        # block coordinates are zero-based, left inclusive, and right-end exclusive
        alignment_segments[0][0] += 1 # adjust for zero-based.  note, end position doesn't need to be adjusted. 
        
        for aligned_pair in aligned_pairs:
            aligned_pair = list(aligned_pair)
            aligned_pair[0] += 1
            
            # extend earlier stored segment or append new one
            if aligned_pair[0] - alignment_segments[-1][1] < self._read_aln_gap_merge_int:
                # extend rather than append
                alignment_segments[-1][1] = aligned_pair[1]
            else:
                # append, as too far apart from prev
                alignment_segments.append(list(aligned_pair))


        # trim short terminal segments from each end
        while (len(alignment_segments) > 1 and
            alignment_segments[0][1] - alignment_segments[0][0] + 1 < self._min_terminal_splice_exon_anchor_length):

            alignment_segments.pop(0)

        while (len(alignment_segments) > 1 and
            alignment_segments[len(alignment_segments)-1][1] - alignment_segments[len(alignment_segments)-1][0] + 1 < self._min_terminal_splice_exon_anchor_length):

            alignment_segments.pop()

            
        return alignment_segments
    
