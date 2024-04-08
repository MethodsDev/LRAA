

test:
	../../LRAA --genome data/SIRVs1-7.genome.fa --bam data/sim.fasta.mm2.bam --norm 0
	gffcompare -r data/SIRVs1-7.annot.gtf LRAA.gtf -o gffcompare
	./eval/eval_sirvs.py data/SIRVs1-7.gene_trans_map gffcompare.LRAA.gtf.refmap
	diff gffcompare.LRAA.gtf.refmap.full.summary.tsv data/latest_summary.tsv



