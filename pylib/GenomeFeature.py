#!/usr/bin/env python
# encoding: utf-8

import sys, os, re
import logging
import string
import pysam
from collections import defaultdict


logger = logging.getLogger(__name__)


class GenomeFeature:

    def __init__(self, contig_acc, lend, rend):
        self._contig_acc = contig_acc
        self._lend = lend
        self._rend = rend
        self._id = "__id_not_set__"

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
    
    def get_id(self):
        return self._id

        
class Intron(GenomeFeature):

    intron_id_counter = 0
    
    def __init__(self, contig_acc, lend, rend, orient, count):
        super().__init__(contig_acc, lend, rend)
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

    def __init__(self, contig_acc, lend, rend, mean_coverage):
        super().__init__(contig_acc, lend, rend)
        self._mean_coverage = mean_coverage

        Exon.exon_id_counter += 1
        self._id = "E:{}".format(Exon.exon_id_counter)

        return

    def get_read_support(self):
        return self._mean_coverage
    
    def __repr__(self):
        return("Exon: {} {}-{} mean_cov:{}".format(self._id, self._lend, self._rend, self._mean_coverage))


