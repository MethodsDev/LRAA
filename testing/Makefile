

run_small:
	../pasa-salraa --bam alignments.1-50k.bam --genome minigenome.fa --CPU 4 --outputfile pasa-salraa.small.gtf


run_large:
	../pasa-salraa --bam alignments.bam --genome minigenome.fa --CPU 4 --outputfile pasa-salraa.large.gtf


run_illum_n_pbio:
	cd illum_n_pbio/ && ../../pasa-salraa  --bam illum_n_pacbio.bam --genome illum_n_pacbio.genome.fa --CPU 4


clean:
	./cleanme.pl

