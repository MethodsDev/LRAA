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

logger = logging.getLogger(__name__)

class PASA_SALRAA:


    def __init__(self, splice_graph):

        self._splice_graph = splice_graph

        return


    def populate_read_multi_paths(self, contig_acc, contig_seq, bam_file):

         bam_extractor = Bam_alignment_extractor(bam_file)
         pretty_alignments = bam_extractor.get_read_alignments(contig_acc, pretty=True)

         grouped_alignments = self._group_alignments_by_read_name(pretty_alignments)

         print(grouped_alignments)

         return

         

    def _group_alignments_by_read_name(self, pretty_alignments):

        grouped_alignments = defaultdict(list)

        for pretty_alignment in pretty_alignments:
            pysam_alignment = pretty_alignment.get_pysam_alignment()
            read_name = pysam_alignment.query_name
            grouped_alignments[read_name].append(pretty_alignment)

        return grouped_alignments
    
