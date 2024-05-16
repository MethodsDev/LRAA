#!/usr/bin/env python
# encoding: utf-8

import sys, os, re
import subprocess
import logging
import string
import pysam
from collections import defaultdict


logger = logging.getLogger(__name__)


class Pretty_alignment:

    def __init__(self, pysam_alignment, pretty_alignment_segments):

        self._pysam_alignment = pysam_alignment
        self._pretty_alignment_segments = pretty_alignment_segments # defined in Bam_alignment_extractor //TODO: move logic here.

        #if pysam_alignment.has_tag("RG") and pysam_alignment.get_tag("RG") == "PBLR":
        #    self._read_type = "PBLR"
        #else:
        #    self._read_type = "ILMN"
        self._read_type = "PacBio"

        self.left_soft_clipping, self.right_soft_clipping = self._get_read_soft_clipping_info(pysam_alignment)

        
        
    def __repr__(self):
        return str(self._pretty_alignment_segments)

    def get_read_name(self):
        return self._pysam_alignment.query_name
    
        
    def get_pysam_alignment(self):
        return self._pysam_alignment

    def get_pretty_alignment_segments(self):
        return self._pretty_alignment_segments


    def get_alignment_span(self):
        lend = self._pretty_alignment_segments[0][0]
        rend = self._pretty_alignment_segments[-1][1]

        return(lend, rend)
    

    def get_read_type(self):
        return self._read_type

    
    def _get_read_soft_clipping_info(self, pysam_alignment=None):

        if pysam_alignment is None:
            pysam_alignment = self._pysam_alignment
        
        cigar_tuples = pysam_alignment.cigartuples

        S=4 # soft clipping cigar code in pysam
        
        left_soft_clipping = cigar_tuples[0][1] if cigar_tuples[0][0] == S else 0

        right_soft_clipping = cigar_tuples[-1][1] if cigar_tuples[-1][0] == S else 0

        return left_soft_clipping, right_soft_clipping


    

    def has_soft_clipping(self):
        left_soft_clip, right_soft_clip = self._get_read_soft_clipping_info()
        if left_soft_clip > 0 or right_soft_clip > 0:
            return True
        else:
            return False



    @classmethod
    def try_correct_alignments(cls, pretty_alignments_list, splice_graph, contig_seq):
    

        for pretty_alignment in pretty_alignments_list:

            if not pretty_alignment.has_soft_clipping():
                continue
            
            alignment_segments = pretty_alignment.get_pretty_alignment_segments()

            left_soft_clipping, right_soft_clipping = pretty_alignment._get_read_soft_clipping_info()

            read_sequence = pretty_alignment._pysam_alignment.query_sequence
            
            if left_soft_clipping:
                left_alignment_segment = alignment_segments[0]
                exon_seg_lend, exon_seg_rend = left_alignment_segment
                overlapping_introns = list()
                for overlapping_intron in splice_graph.get_overlapping_introns(exon_seg_lend, exon_seg_rend):
                    intron_lend, intron_rend = overlapping_intron.get_coords()
                    if intron_rend > exon_seg_lend and intron_rend < exon_seg_rend:
                        overlapping_introns.append(overlapping_intron)

                print("Got overlapping introns: {}".format(overlapping_introns))


            if right_soft_clipping:
                right_alignment_segment = alignment_segments[-1]



        return

    
        
