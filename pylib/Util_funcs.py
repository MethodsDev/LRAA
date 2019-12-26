#!/usr/bin/env python
# encoding: utf-8

import sys, os, re
import subprocess
import logging
import string

from collections import defaultdict

logger = logging.getLogger(__name__)


def retrieve_contig_seq_from_fasta_file(fasta_filename, contig_acc):

    contig_seq_str = subprocess.check_output("samtools faidx {} {}".format(contig_acc, fasta_filename),
                                             shell=True,
                                             encoding="utf-8")
    
    contig_seq_str = contig_seq_str.upper()
    contig_seq_str = contig_seq_str.split("\n")
    contig_seq_str = contig_seq_str[1:]
    contig_seq_str = "".join(contig_seq_str)

    contig_seq_str = re.sub("\s", "", contig_seq_str) # just in case

    return(contig_seq_str)


def coordpairs_overlap(coordset_A, coordset_B):

    A_lend, A_rend = coordset_A

    B_lend, B_rend = coordset_B

    if A_lend <= B_rend and A_rend >= B_lend:
        return True
    else:
        return False

    
