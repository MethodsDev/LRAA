
## global vars / constants

SPACER="???"

DEBUG=False

config = {
    
    # read alignment criteria
    'min_per_id' : 98,
    'min_mapping_quality' : 20,
    'try_correct_alignments' : False,
    'max_softclip_realign_test' : 10,
    
    # splice graph construction criteria
    'min_alt_splice_freq' : 0.01,
    'min_alt_unspliced_freq' : 0.01,
    'min_feature_frac_overlap' : 0.50,
    'max_exon_spur_length' : 5, # exon spurs not tied to TSS or PolyA and at most this length get pruned
    
    ## TSS config
    'infer_TSS' : True, # include TSS feature in read path assignments
    'max_dist_between_alt_TSS_sites' : 0,
    'min_alignments_define_TSS_site' : 3,
    'max_frac_compatible_expression' : 0.5, # if 2 transcripts (A,B) are compatible and A contains B, and B < this fraction expression, B excluded as likely degradation product
    'min_TSS_iso_fraction' : 0.01, # like min_isoform_fraction, but specifically targeting TSS features
    'max_frac_alt_TSS_from_degradation' : 0.20,
    'max_soft_clip_at_TSS' : 3,
    
    ## polyA site config
    'infer_PolyA' : True, # include PolyA site feature in read path assignments
    'max_dist_between_alt_polyA_sites' : 50,
    'min_alignments_define_polyA_site' : 3,
    


    ## read assignment to transcript criteria
    'fraction_read_align_overlap' : 0.75, # min fraction of read length that must overlap the compatible transcript isoform structure
    
    
    # misc settings
    'min_path_score' : 2, # min number of reads required for reporting isoform

    
    # transcript criteria
    'min_transcript_length' : 200,
    'min_isoform_fraction' : 0.01,
    
    
    # assembly
    'normalize_max_cov_level' : 10000,
    'restrict_asm_to_collapse' : True, # if True, no chaining of overlapping/extended paths
    'collapse_alt_TSS_and_PolyA' : True  # if True, collapses paths that are overlapping and contained but differ in TSS or PolyA

    
    }

