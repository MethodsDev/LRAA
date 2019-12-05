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

        

    
    
