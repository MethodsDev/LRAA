#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;


## we delete all files we don't need in this directory. Be careful in case users try running it somewhere else, outside this dir.
chdir $FindBin::Bin or die "error, cannot cd to $FindBin::Bin";



my @files_to_keep = qw (


alignments.1-50k.bam
alignments.1-50k.bam.bai
alignments.bam
alignments.bam.bai
minigenome.fa
minigenome.fa.fai
minigenome.gtf
runMe.sh
runMe.small.sh
cleanme.pl



);

my %keep = map { + $_ => 1 } @files_to_keep;


foreach my $file (<*>) {
	
	if (-f $file && ! $keep{$file}) {
		print STDERR "-removing file: $file\n";
		unlink($file);
	}
}


exit(0);
