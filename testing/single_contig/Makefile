
test_reads_only:
	../../LRAA --bam pacbio.PBLR.bam --genome ref_genome.fa --CPU 4 --outputfile LRAA.gtf  --normalize_max_cov_level 0


test_include_ref_annot_guide:
	../../LRAA --bam pacbio.PBLR.bam --genome ref_genome.fa --CPU 4 --outputfile LRAA.guided.gtf  --normalize_max_cov_level 0 --gtf ref_annot.gtf


all: test_reads_only test_include_ref_annot_guide


clean:
	rm -rf ./__* ./LRAA*.gtf


