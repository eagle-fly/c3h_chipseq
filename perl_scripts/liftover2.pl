
#!/usr/bin/perl -w
use strict;

# convert the track_data file from mm9 to mm10 genome
# Author: Ying Zhang; yzhang@bccrc.ca
# Date: Jan 14, 2015

my $old_track = '../documents/ListIAPinC3H&NOTB6mm9.txt';
my $old_loci = '../documents/ListIAPinC3H&NOTB6_loci_mm9.txt';
my $new_loci = '../documents/ListIAPinC3H&NOTB6mm10.txt';

open (OLDTRACK, "$old_track") or die "cannot open track data file: $!";
open (OLDLOCI, ">$old_loci") or die "cannot open track data file: $!";

my %data;
while (<OLDTRACK>) {
    my ($chr, $start, $end, @other_fields) = split /\t/, $_;
    print OLDLOCI "chr$chr\t$start\t$end\n";
}

close OLDTRACK;
close OLDLOCI;

# convert old loci from mm10 to mm9 by using the liftOver tool
my $liftover = "./liftOver \"$old_loci\" mm9ToMm10.over.chain \"$new_loci\" ../documents/mm10.unmaped.txt";
system ($liftover) and die "cannot liftOver\n";


print "\nDone!\n";