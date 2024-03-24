

run_small:
	../pasa-salraa --bam alignments.1-50k.bam --genome minigenome.fa --CPU 4 --outputfile pasa-salraa.small.gtf


specific_region_only:
	../pasa-salraa --bam alignments.1-50k.bam --genome minigenome.fa --CPU 4 --outputfile pasa-salraa.small.restricted_range.gtf --restrict_to_contig minigenome --restrict_region_range 24920-29235

run_large:
	../pasa-salraa --bam alignments.bam --genome minigenome.fa --CPU 4 --outputfile pasa-salraa.large.gtf



clean:
	./cleanme.pl

