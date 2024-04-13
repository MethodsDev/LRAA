#!/bin/bash

set -ex

gffcompare  -r SIRV5.annot.gtf LRAA.gtf -o gffcompare

../eval/eval_sirvs.py SIRV5.gene_trans_map gffcompare.LRAA.gtf.refmap
