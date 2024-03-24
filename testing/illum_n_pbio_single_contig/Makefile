

pacbio:
	../../pasa-salraa --bam pacbio.PBLR.bam --genome ref_genome.fa --CPU 4 --outputfile pasa-salraa.pbio.gtf --mpm_monitor


illumina:
	../../pasa-salraa --bam illumina.bam --genome ref_genome.fa --CPU 4 --outputfile pasa-salraa.illumina.gtf


both:
	../../pasa-salraa  --bam illum_n_pacbio.bam --genome ref_genome.fa --CPU 4



clean:
	rm -rf ./__* ./pasa-salraa*
