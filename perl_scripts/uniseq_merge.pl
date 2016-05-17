
#!/usr/bin/perl -w
use strict;

# Consolidate the two uniseq files for both read ends into a single file
# by removing the shorter uniseq when two mate reads are present with uniseq
# Author: Ying Zhang; yzhang@bccrc.ca
# Date: Nov 13, 2014

my $dataset = "c3h.input";
my $uniseq_file1 = "../results/uniseq/$dataset.1.uniseq";
my $uniseq_file2 = "../results/uniseq/$dataset.2.uniseq";

open (UNISEQ1, "$uniseq_file1") or die "cannot open uniseq file: $!";
my %uniseq;
while (<UNISEQ1>) {
    chomp;
    my ($id, $seq) = split /\t/;
    $uniseq{$id} = $seq;
}
open (UNISEQ2, "$uniseq_file2") or die "cannot open uniseq file: $!";
while (<UNISEQ2>) {
    chomp;
    my ($id, $seq) = split /\t/;
    my $mate_id = $id;
    substr $mate_id, -1, 1, '1';
    if (exists $uniseq{$mate_id}) {
        if (length $uniseq{$mate_id} < length $seq) {
            $uniseq{$id} = $seq;
            delete $uniseq{$mate_id};
        }
    } else {
        $uniseq{$id} = $seq;
    }    
}

my $output_file = "../results/uniseq/$dataset.uniseq";
open (OUTFILE, ">$output_file") or die "cannot open OUTPUT file: $!";
foreach my $id (sort keys %uniseq) {
    print OUTFILE "$id\t$uniseq{$id}\n";
}
print "Done!\n";
