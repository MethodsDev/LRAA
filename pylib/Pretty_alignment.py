#!/usr/bin/env pythonOA
# encoding: utf-8

import sys, os, re
import subprocess
import logging
import string
import pysam
import PASA_SALRAA_Globals
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

        logger.info("Attempting to correct alignments at soft-clips")
        
        max_softclip_realign_test = PASA_SALRAA_Globals.config['max_softclip_realign_test']
        
        
        for pretty_alignment in pretty_alignments_list:

            if not pretty_alignment.has_soft_clipping():
                continue
            
            alignment_segments = pretty_alignment.get_pretty_alignment_segments()

            left_soft_clipping, right_soft_clipping = pretty_alignment._get_read_soft_clipping_info()

            read = pretty_alignment._pysam_alignment
            
            read_sequence = read.query_sequence

            # get mapping of genome pos -> read pos
            aligned_pairs = dict([ (y+1, x+1) for x,y in read.get_aligned_pairs(matches_only=True) ])


            # ignore soft-clips involving just polyA - only useful if sim data here with tacked-on polyA - could make more useful, but just trim polyA before running is easiest.
            if (read.is_forward and
                right_soft_clipping > 0 and
                right_soft_clipping <= max_softclip_realign_test and
                re.match("^A+$", read_sequence[ (-1*right_soft_clipping):], re.I) is not None):
                continue

            
            if (read.is_reverse and
                left_soft_clipping > 0 and
                left_soft_clipping <= max_softclip_realign_test and
                re.match("^T+$", read_sequence[0:left_soft_clipping], re.I) is not None):
                continue

            
            if left_soft_clipping > 0 and left_soft_clipping <= max_softclip_realign_test:
                left_alignment_segment = alignment_segments[0]
                exon_seg_lend, exon_seg_rend = left_alignment_segment
                overlapping_introns = list()
                for overlapping_intron in splice_graph.get_overlapping_introns(exon_seg_lend, exon_seg_rend):
                    intron_lend, intron_rend = overlapping_intron.get_coords()
                    if intron_rend > exon_seg_lend and intron_rend < exon_seg_rend:
                        overlapping_introns.append(overlapping_intron)

                print("Got left overlapping introns: {}".format(overlapping_introns))
                for left_overlapping_intron in overlapping_introns:
                    intron_lend, intron_rend = left_overlapping_intron.get_coords()

                    intron_adjacent_pos = intron_rend + 1
                    if intron_adjacent_pos not in aligned_pairs:
                        continue
                    read_rend = aligned_pairs[intron_adjacent_pos] - 1
                    if read_rend -1 <=  max_softclip_realign_test:
                        left_read_seq = read_sequence[0:read_rend]
                        print("Checking read realignment for {}".format(left_read_seq))
                        
                        genomic_rend = intron_lend - 1
                        genomic_lend = genomic_rend - len(left_read_seq) + 1
                        
                        genomic_substr = contig_seq[genomic_lend - 1 : genomic_rend]
                        print("Comparing to genomic rend seq: {}".format(genomic_substr))
                        
                        if left_read_seq.upper() == genomic_substr.upper():
                            print("\tLeft MATCH FOUND")
                            # do reassignment:
                            alignment_segments[0][0] = intron_rend + 1
                            alignment_segments.insert(0, [genomic_lend, genomic_rend])
                            break
            
            
            if right_soft_clipping > 0 and right_soft_clipping <= max_softclip_realign_test:
                right_alignment_segment = alignment_segments[-1]
                exon_seg_lend, exon_seg_rend = right_alignment_segment
                overlapping_introns = list()
                for overlapping_intron in splice_graph.get_overlapping_introns(exon_seg_lend, exon_seg_rend):
                    intron_lend, intron_rend = overlapping_intron.get_coords()
                    if intron_lend > exon_seg_lend and intron_lend < exon_seg_rend:
                        overlapping_introns.append(overlapping_intron)

                print("Got right overlapping introns: {}".format(overlapping_introns))
                for right_overlapping_intron in overlapping_introns:
                    intron_lend, intron_rend = right_overlapping_intron.get_coords()
                    intron_adjacent_pos = intron_lend - 1
                    if intron_adjacent_pos not in aligned_pairs:
                        continue
                    read_lend = aligned_pairs[intron_adjacent_pos]
                    if len(read_sequence) - read_lend + 1 <= max_softclip_realign_test:
                        right_read_seq = read_sequence[read_lend:]
                        print("Checking read realignment for {}".format(right_read_seq))
                        genomic_lend = intron_rend + 1
                        genomic_rend = genomic_lend + len(right_read_seq) - 1
                        genomic_substr = contig_seq[genomic_lend - 1 : genomic_rend]
                        print("Comparing to genomic rend seq: {}".format(genomic_substr))
                        if right_read_seq.upper() == genomic_substr.upper():
                            print("\tRight MATCH FOUND")
                            # do reassignment:
                            alignment_segments[-1][1] = intron_lend -1
                            alignment_segments.append([genomic_lend, genomic_rend])
                            break

        return

    
        
