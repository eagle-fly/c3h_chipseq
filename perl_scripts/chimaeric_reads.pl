
#!/usr/bin/perl -w
use strict;

# get the full read sequences of the selected chimaeric reads in mappable_uniseq.fa from the masked reads
# Author: Ying Zhang; yzhang@bccrc.ca
# Date: Nov 13, 2014

my $dataset = "c3h.input";
my $mappable_file = "../results/uniseq/$dataset.mappable.fa";
my $masked_file1 = "../data/rm_out/$dataset.1/$dataset.1.fa.masked";
my $masked_file2 = "../data/rm_out/$dataset.2/$dataset.2.fa.masked";

open (MAPPABLE, "$mappable_file") or die "cannot open mappable file: $!";
open (MASKED1, "$masked_file1") or die "cannot open masked file: $!";
open (MASKED2, "$masked_file2") or die "cannot open masked file: $!";

my $output_file = "../results/uniseq/chimaeric_reads.fa";
open (OUTFILE, ">$output_file") or die "cannot open OUTPUT file: $!";

while (my $seq_id = <MAPPABLE>) {
    chomp $seq_id;
    <MAPPABLE>; # skip the sequence line
    
    my $m = substr $seq_id, -1;
    if ($m == 1) {
        while (my $read_id = <MASKED1>) {
            chomp $read_id;
            my $read_seq = <MASKED1>;
            chomp $read_seq;
            $read_seq .= <MASKED1>;
            chomp $read_seq;
            $read_seq = uc $read_seq;
            if ($seq_id eq $read_id) {
                print OUTFILE "$read_id\n";
                print OUTFILE "$read_seq\n";
                last;
            }
        }
    } else {
        while (my $read_id = <MASKED2>) {
            chomp $read_id;
            my $read_seq = <MASKED2>;
            chomp $read_seq;
            $read_seq .= <MASKED2>;
            chomp $read_seq;
            $read_seq = uc $read_seq;
            if ($seq_id eq $read_id) {
                print OUTFILE "$read_id\n";
                print OUTFILE "$read_seq\n";
                last;
            }
        }
    }
}
print "Done!\n";

