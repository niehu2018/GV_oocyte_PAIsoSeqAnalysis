#!/usr/bin/perl -w

use strict;

my $usage = "perl $0 ccs.bam > ccs.pass.txt\n";
if ( @ARGV != 1 ) {
    print $usage;
    exit;
}

my $infile = $ARGV[0];
open IN, "samtools view $infile | " or die "Could not open file $infile $!\n";
while (<IN>) {
    chomp;
    my @fields = split;
    my $id     = $fields[0];
    my $pass   = 0;
    if (/\s+np:i:(\d+)\s+/) {
        $pass = $1;
    }
    print "$id\t$pass\n" if $pass > 0;
}
close IN;
