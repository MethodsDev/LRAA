#!/usr/bin/env python
# encoding: utf-8

import sys, os, re
import subprocess
import logging
import string
import pysam

logger = logging.getLogger(__name__)


_MAX_CONTIG_LENGTH = 1e10
_read_aln_gap_merge_int = 10


class splice_graph:

    _genome_fasta_filename = ""
    _alignments_bam_filename = ""

    _contig_acc = ""
    _contig_seq_str = ""
    _contig_base_cov = []
    
    _splice_graph = dict()

    
    
    
    def __init__(self):
        
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
        if (contig_len > _MAX_CONTIG_LENGTH):
            raise RuntimeError("contig length {} exceeds maximum allowed {}".format(contig_len, _MAX_CONTIG_LENGTH))

        
        # init depth of coverage array
        self._contig_base_cov = [0 for i in range(0,contig_len)]
        
        samfile = pysam.AlignmentFile(alignments_bam_file, "rb")
        for read in samfile.fetch(contig_acc):
            print(read)
            aligned_pairs = read.get_blocks()
            print(aligned_pairs)


            
        
        return

    

    def _retrieve_contig_seq(self):

        contig_seq_str = subprocess.check_output("samtools faidx {} {}".format(self._genome_fasta_filename,
                                                                               self._contig_acc),
                                                 shell=True,
                                                 encoding="utf-8")
        
        contig_seq_str = contig_seq_str.upper()
        contig_seq_str = re.sub("\s", "", contig_seq_str)
        
        self._contig_seq_str = contig_seq_str

        # print(contig_seq_str)
        
        return(contig_seq_str)



    
        
