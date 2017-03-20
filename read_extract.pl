#!/usr/bin/env perl
use strict;
use warnings;

# TODO: beatify the following comments.
# Custome document for extracting specific reads (don't have N at specific 
# location) from fastq
# Usage: perl seq_extract.pl [fastqfile] [start] [length] [number of reads]

my $file = $ARGV[0];
my $fastq;
if ($file =~ /.gz$/) {
    open($fastq, "-|", "gzip -dc $file") or die "Open gzip file failed: $!"
}
else {
    open($fastq,  "<".$file)  or die "Can't open $file: $!";
}

my ($seq_id, $seq, $subseq, $pqs_id, $pqs);
my $count = 0;

while (<$fastq>) { # assigns each line in turn to $_
	$seq_id = $_;
	$seq = <$fastq>;
	chomp $seq;
	$subseq = substr $seq, $ARGV[1], $ARGV[2];
	$pqs_id = <$fastq>;
	$pqs = <$fastq>;
    if ($subseq =~ m/[nN]/) {
        print $seq_id;
        print $seq."\n";
        print $pqs_id;
        print $pqs;
        $count++
    }
    if ($count > $ARGV[3] - 1) {
        last;
    }
}

undef $fastq;
