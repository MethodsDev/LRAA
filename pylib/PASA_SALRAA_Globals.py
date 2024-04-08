
## global vars / constants

SPACER="???"

DEBUG=False

config = {
    
    # read alignment criteria
    'min_per_id' : 98,


    # splice graph construction criteria

    'infer_TSS' : True, # include TSS feature in read path assignments
    'infer_PolyA' : True, # include PolyA site feature in read path assignments
    
    'max_dist_between_alt_TSS_sites' : 50,
    'max_dist_between_alt_polyA_sites' : 50,

    'min_alignments_define_TSS_site' : 3,
    'min_alignments_define_polyA_site' : 3,

    'max_frac_alt_TSS_from_degradation' : 0.20,
    
    # misc settings
    'min_path_score' : 1,

    
    
    

    # transcript criteria
    'min_transcript_length' : 200,

    'min_isoform_fraction' : 0.01
        
    }

