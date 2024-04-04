#!/usr/bin/env python3

import sys, os, re
import Transcript
import MultiPath
import MultiPathCounter
import Simple_path_utils as SPU
from collections import defaultdict
from PASA_SALRAA_Globals import SPACER, DEBUG
import logging


logger = logging.getLogger(__name__)



class Quantify:

    def __init__(self):

        self._path_node_id_to_gene_ids = defaultdict(set)

        self._gene_id_to_transcript_objs = defaultdict(set)

        self._read_name_to_multipath = dict()
        
        return


    def quantify(self, splice_graph, transcripts, mp_counter, run_EM=True):

        assert type(transcripts) == list
        assert type(transcripts[0]) == Transcript.Transcript
        assert type(mp_counter) == MultiPathCounter.MultiPathCounter

        
        # assign path nodes to gene
        # also assign gene_id to transcript objs
        self._assign_path_nodes_to_gene(transcripts)

        self._assign_reads_to_transcripts(splice_graph, mp_counter)

        transcript_to_read_count = self._estimate_isoform_read_support(transcripts, run_EM)

        return transcript_to_read_count

        
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



    def _assign_reads_to_transcripts(self, splice_graph, mp_counter, fraction_read_align_overlap=0.75):

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
            transcripts_assigned = self._assign_path_to_transcript(splice_graph, mp, gene_isoforms, fraction_read_align_overlap)
            if transcripts_assigned is None:
                logger.debug("mp_count_pair {} maps to gene but no isoform(transcript)".format(mp_count_pair))
            else:
                logger.debug("mp_count_pair {} maps to transcripts: {}".format(mp_count_pair, transcripts_assigned))
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
        

        
    def _assign_path_to_transcript(self, splice_graph, mp, transcripts, fraction_read_align_overlap):
        
        assert type(mp) == MultiPath.MultiPath
        assert type(transcripts) == set, "Error, type(transcripts) is {} not set ".format(type(transcripts))
        assert type(list(transcripts)[0]) == Transcript.Transcript
        assert fraction_read_align_overlap >= 0 and fraction_read_align_overlap <= 1.0, "Error, fraction_read_align_overlap must be between 0 and 1.0"
        
        read_sp = mp.get_simple_path()

        # store read name to mp for later debugging.
        for read_name in mp.get_read_names():
            self._read_name_to_multipath[read_name] = mp

        
        transcripts_compatible_with_read = list()
        
        for transcript in transcripts:
            transcript_sp = transcript._simplepath
            assert transcript_sp is not None

            if (SPU.are_overlapping_and_compatible_NO_gaps_in_overlap(transcript_sp, read_sp)
                and
                SPU.fraction_read_overlap(splice_graph, read_sp, transcript_sp) >= fraction_read_align_overlap):
                
                #print("Read {} compatible with transcript {}".format(read_sp, transcript_sp))
                transcripts_compatible_with_read.append(transcript)
                
                        
        return transcripts_compatible_with_read



    def _estimate_isoform_read_support(self, transcripts, run_EM):


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
            logger.debug(f"-assigning transcript {transcript_id} read count: {transcript_read_count_total}")

            
        ## DEBUGGING
        logger.debug("# Isoform read assignments:\n")
        for read_name in read_name_to_transcripts:
            transcripts_read_assigned = read_name_to_transcripts[read_name]
            logger.debug("read_name {} assigned to {}".format(read_name, transcripts_read_assigned))
            if len(transcripts_read_assigned) > 1:
                logger.debug("*** Splitting read: {} across {} transcripts: {}".format(read_name, len(transcripts_read_assigned), transcripts_read_assigned))
            
     
        if run_EM:

            ## go through multiple rounds of EM

            for i in range(1, 100):
                
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
            

        return transcript_to_read_count
        
    
                
    def report_quant_results(self, transcripts, transcript_to_read_count, ofh_quant_vals, ofh_read_tracking):
                
        ## generate final report.
        for transcript in transcripts:
            transcript_id = transcript.get_transcript_id()
            gene_id = transcript.get_gene_id()
            counts = transcript_to_read_count[transcript_id]

            readnames = transcript.get_read_names()
            readnames = sorted(readnames)

            logger.info("\t".join([gene_id, transcript_id, f"{counts:.1f}"]))
            print("\t".join([gene_id, transcript_id, f"{counts:.1f}"]), file=ofh_quant_vals)

            if (DEBUG):
                print("transcript_id\t{}\n{}".format(transcript_id, transcript._simplepath), file=ofh_read_tracking)
                for readname in readnames:
                    print("read:\t{}\n{}".format(readname, self._read_name_to_multipath[readname].get_simple_path()), file=ofh_read_tracking)
                print("\n", file=ofh_read_tracking)
            else:
                print("\t".join([gene_id, transcript_id, ",".join(readnames)]), file=ofh_read_tracking)
            
        return 
    


    @staticmethod
    def filter_isoforms_by_min_isoform_fraction(transcripts, min_isoform_fraction, run_EM):

        logger.info("Filtering transcripts according to min isoform fraction: {}".format(min_isoform_fraction))
        
        isoforms_were_filtered = True # init for loop

        q = Quantify()


        filtering_round = 0
        
        while isoforms_were_filtered:

            filtering_round += 1

            num_filtered_isoforms = 0
            num_total_isoforms = len(transcripts)
            
            transcripts_retained = list()
            
            isoforms_were_filtered = False # update to True if we do filter an isoform out.
            
            # run (or rerun) quant
            transcript_to_read_count = q._estimate_isoform_read_support(transcripts, run_EM)


            gene_to_transcripts = defaultdict(list)
            for transcript in transcripts:
                gene_id = transcript.get_gene_id()
                gene_to_transcripts[gene_id].append(transcript)

            for gene_id in gene_to_transcripts:
                transcripts_of_gene = gene_to_transcripts[gene_id]
                if len(transcripts_of_gene) == 1:
                    transcripts_retained.extend(transcripts_of_gene)
                    continue

                # evaluate isoform fraction
                sum_gene_reads = 0
                for transcript_of_gene in transcripts_of_gene:
                    transcript_read_count = transcript_to_read_count[ transcript_of_gene.get_transcript_id() ]
                    sum_gene_reads += transcript_read_count

                logger.info("gene_id {} has total reads: {}".format(gene_id, sum_gene_reads))
                
                    
                for transcript_of_gene in transcripts_of_gene:
                    transcript_id =  transcript_of_gene.get_transcript_id()
                    transcript_read_count = transcript_to_read_count[ transcript_id ]
                    isoform_frac = transcript_read_count / sum_gene_reads
                    logger.info("\ttranscript_id {} has {} reads = {} isoform fraction of {}".format(
                        transcript_id,
                        transcript_read_count,
                        isoform_frac,
                        gene_id))
                    
                    if isoform_frac < min_isoform_fraction:                        
                        isoforms_were_filtered = True
                        num_filtered_isoforms += 1
                    else:
                        transcripts_retained.append(transcript_of_gene)


            logger.info("isoform filtering round {} involved filtering of {} isoforms / {} total isoforms of {} genes".format(
                filtering_round,
                num_filtered_isoforms,
                num_total_isoforms,
                len(gene_to_transcripts)) )
            
            # reset list of transcripts
            transcripts = transcripts_retained


            
        return transcripts





        
