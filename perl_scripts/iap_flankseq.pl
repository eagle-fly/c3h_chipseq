
#!/usr/bin/perl -w
use strict;

# Get the flanking unique sequence from reads containing IAP
# Author: Ying Zhang; yzhang@bccrc.ca
# Date: Nov 6, 2014

my $dataset = "c3h.input.2";
my $masked_file = "../data/rm_out/$dataset/$dataset.fa.masked";
my $iap_file = "../data/rm_out/$dataset/$dataset.iap";
my $output_file = "../results/uniseq/$dataset.uniseq";

open (MASKED, "$masked_file") or die "cannot open masked file: $!";
open (IAPFILE, "$iap_file") or die "cannot open IAP file: $!";
open (OUTFILE, ">$output_file") or die "cannot open OUTPUT file: $!";

while (<IAPFILE>) {
    chomp;
    my @fields = split /\s+/;
    if ($fields[10] =~ "-int") {    # ignore IAP internal sequences
        next;
    }
    my $read_id = $fields[5];
    
    # find the corresponding read in the masked sequence file
    while (my $seq_id = <MASKED>) {
        chomp $seq_id;
        $seq_id = substr $seq_id, 1;
        my $seq = <MASKED>;
        chomp $seq;
        $seq .= <MASKED>;
        chomp $seq;
        
        my $uni_seq;
        if ($seq_id eq $read_id) {
            $uni_seq = longest_ucstr($seq);
            print OUTFILE "$read_id\t$uni_seq\n";
            last; # stop checking the MASKED file
        }
    }
    
}
close IAPFILE;
close MASKED;
print "Done!\n";

sub longest_ucstr {
    my $seq = shift;
    
    my @uc_all = $seq =~ /[ACGT]+/g;
    my $longest = "";
    foreach my $uc (@uc_all) {
        if (length($longest) < length($uc)) {
            $longest = $uc;
        }        
    }
    return $longest;
}