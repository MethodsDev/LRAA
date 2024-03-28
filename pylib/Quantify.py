#!/usr/bin/env python3

import sys, os, re
import Transcript
import MultiPath
import MultiPathCounter
import Simple_path_utils as SPU
from collections import defaultdict
from PASA_SALRAA_Globals import SPACER
import logging


logger = logging.getLogger(__name__)



class Quantify:

    def __init__(self):

        self._path_node_id_to_gene_ids = defaultdict(set)

        self._gene_id_to_transcript_objs = defaultdict(set)
        
        return


    def quantify(self, transcripts, mp_counter):

        assert type(transcripts) == list
        assert type(transcripts[0]) == Transcript.Transcript
        assert type(mp_counter) == MultiPathCounter.MultiPathCounter

        
        # assign path nodes to gene
        # also assign gene_id to transcript objs
        self._assign_path_nodes_to_gene(transcripts)

        self._assign_reads_to_transcripts(mp_counter)


    def _assign_path_nodes_to_gene(self, transcripts):

        for transcript in transcripts:
            
            simplepath = transcript._simplepath
            assert simplepath is not None, "Error, simplepath not set for transcript obj"

            transcript_id = transcript.get_transcript_id()
            gene_id = transcript.get_gene_id()
            self._gene_id_to_transcript_objs[gene_id].add(transcript)
            
            for node_id in simplepath:
                if node_id != SPACER:
                    self._path_node_id_to_gene_ids[node_id].add(gene_id)

        return



    def _assign_reads_to_transcripts(self, mp_counter):

        # assign to gene based on majority voting of nodes.
        # TODO:// might want or need this to involve length and/or feature type weighted shared node voting

        mp_count_pairs = mp_counter.get_all_MultiPathCountPairs()

        gene_unanchored_mp_count_pairs = list()


        num_paths_total = 0
        num_read_counts_total = 0

        num_paths_anchored_to_gene = 0
        num_read_counts_anchored_to_gene = 0
        
        num_paths_assigned = 0
        num_read_counts_assigned = 0

        
        
        for mp_count_pair in mp_count_pairs:
            mp, count = mp_count_pair.get_multipath_and_count()

            num_paths_total += 1
            num_read_counts_total += count
            
            sp = mp.get_simple_path()

            top_gene = self._get_gene_with_best_node_matches_to_simplepath(sp)

            if top_gene is None:
                gene_unanchored_mp_count_pairs.append(mp_count_pair)
                logger.debug("mp_count_pair unanchored: " + str(mp_count_pair))
                
                continue

            logger.debug("mp_count_pair {} anchored to gene: {}".format(mp_count_pair, top_gene))
            num_paths_anchored_to_gene += 1
            num_read_counts_anchored_to_gene += count
            
            
            ## assign reads to transcripts
            gene_isoforms = self._gene_id_to_transcript_objs[top_gene]
            transcripts_assigned = self._assign_path_to_transcript(mp, gene_isoforms)
            if transcripts_assigned is None:
                logger.debug("mp_count_pair {} maps to gene but no isoform(transcript)".format(mp_count_pair))
            else:
                for transcript in transcripts_assigned:
                    transcript.add_read_names(mp.get_read_names())
                    
                num_paths_assigned += 1
                num_read_counts_assigned += count

        ## audit summary
        audit_txt = "\n".join(
            [
             "num_paths_total: {}, num_read_counts_total: {}".format(num_paths_total, num_read_counts_total),
             "\tnum_paths_anchored_to_gene: {} = {:.2f}%, num_read_counts_anchored_to_gene: {} = {:.2f}%\n".format(num_paths_anchored_to_gene,
                                                                                                                 num_paths_anchored_to_gene/num_paths_total*100,
                                                                                                                 num_read_counts_anchored_to_gene,
                                                                                                                 num_read_counts_anchored_to_gene/num_read_counts_total*100),
            "\tnum_paths_assigned_to_trans: {} = {:.2f}%, num_read_counts_assigned_to_trans: {} = {:.2f}%\n".format(num_paths_assigned,
                                                                                                                  num_paths_assigned/num_paths_total*100,
                                                                                                                  num_read_counts_assigned,
                                                                                                                  num_read_counts_assigned/num_read_counts_total*100)])
        
                
        logger.info(audit_txt)
                
        return


    def _get_gene_with_best_node_matches_to_simplepath(self, simplepath):

        gene_ranker = defaultdict(int)

        for node in simplepath:
            if node != SPACER:
                if node in self._path_node_id_to_gene_ids:
                    gene_set = self._path_node_id_to_gene_ids[node]
                    for gene_id in gene_set:
                        gene_ranker[gene_id] += 1

        if len(gene_ranker) == 0:
            return None
        else:
            genes_ranked = sorted(gene_ranker.keys(), key=lambda x: gene_ranker[x], reverse=True)
            return genes_ranked[0]
        

        
    def _assign_path_to_transcript(self, mp, transcripts):
        
        assert type(mp) == MultiPath.MultiPath
        assert type(transcripts) == set, "Error, type(transcripts) is {} not set ".format(type(transcripts))
        assert type(list(transcripts)[0]) == Transcript.Transcript

        read_sp = mp.get_simple_path()

        transcripts_compatible_with_read = list()
        
        for transcript in transcripts:
            transcript_sp = transcript._simplepath
            assert transcript_sp is not None

            if SPU.are_overlapping_and_compatible_NO_gaps_in_overlap(transcript_sp, read_sp):
                print("Read {} compatible with transcript {}".format(read_sp, transcript_sp))
                transcripts_compatible_with_read.append(transcript)

        return transcripts_compatible_with_read

    
                
    def report_quant_results(self, transcripts, ofh):


        read_name_to_transcripts = defaultdict(set)
        
        for transcript in transcripts:
            read_names = transcript.get_read_names()
            for read_name in read_names:
                read_name_to_transcripts[read_name].add(transcript)

                
        num_mapped_reads = len(read_name_to_transcripts)


        transcript_to_read_count = defaultdict(float)

        transcript_to_expr_val = defaultdict(float)
        
        ## first round of EM for now - split evenly across mapped transcripts.
        for transcript in transcripts:

            transcript_id = transcript.get_transcript_id()
            
            transcript_read_count_total = 0
            read_names = transcript.get_read_names()
            for read_name in read_names:
                num_transcripts_with_assigned_read = len(read_name_to_transcripts[read_name])
                transcript_read_count_total += 1 / num_transcripts_with_assigned_read
                
            transcript_to_read_count[transcript_id] = transcript_read_count_total
            transcript_to_expr_val[transcript_id] = transcript_read_count_total / num_mapped_reads * 1e6



            
        ## go through multiple rounds of EM
        
        for i in range(1,100):

            logger.info("EM round {}".format(i))
            
            ## fractionally assign reads based on expr values
            transcript_to_read_count.clear()
            
            for transcript in transcripts:
                transcript_id = transcript.get_transcript_id()
                transcript_read_count_total = 0
                read_names = transcript.get_read_names()
                transcript_expr = transcript_to_expr_val[transcript_id] 
                for read_name in read_names:
                    transcripts_with_read = read_name_to_transcripts[read_name]
                    sum_expr = 0
                    for tran_with_read in transcripts_with_read:
                        tran_with_read_id = tran_with_read.get_transcript_id()
                        sum_expr += transcript_to_expr_val[tran_with_read_id]
                    frac_read_assignment = transcript_expr / sum_expr
                    transcript_to_read_count[transcript_id] += frac_read_assignment

            ## recompute expr_vals
            transcript_to_expr_val.clear()

            for transcript in transcripts:
                transcript_id = transcript.get_transcript_id()
                transcript_read_count = transcript_to_read_count[transcript_id]
                transcript_to_expr_val[transcript_id] = transcript_read_count/num_mapped_reads * 1e6
            
                
        ## generate final report.
        for transcript in transcripts:
            transcript_id = transcript.get_transcript_id()
            gene_id = transcript.get_gene_id()
            counts = transcript_to_read_count[transcript_id]
            expr = transcript_to_expr_val[transcript_id]

            print("\t".join([gene_id, transcript_id, f"{counts:.1f}", f"{expr:.3f}"]), file=ofh)


        return 
    
