
#!/usr/bin/perl -w
use strict;

# screen the unique sequence files and the corresponding mate sequence in masked file to evaluate the mappability 
# Author: Ying Zhang; yzhang@bccrc.ca
# Date: Nov 13, 2014

my $dataset = "c3h.input";
my $uniseq_file = "../results/uniseq/$dataset.uniseq";
my $masked_file1 = "../data/rm_out/$dataset.1/$dataset.1.fa.masked";
my $masked_file2 = "../data/rm_out/$dataset.2/$dataset.2.fa.masked";

open (UNISEQ, "$uniseq_file") or die "cannot open uniseq file: $!";
open (MASKED1, "$masked_file1") or die "cannot open masked file: $!";
open (MASKED2, "$masked_file2") or die "cannot open masked file: $!";

my $output_file = "../results/uniseq/$dataset.mappable.fa";
open (OUTFILE, ">$output_file") or die "cannot open OUTPUT file: $!";

while (<UNISEQ>) {
    chomp;
    my ($id, $uniseq) = split /\t/;
    my ($mate_id, $mate_uniseq);
    
    my $m = substr $id, -1;
    if ($m == 1) {
        $mate_id = $id;
        substr $mate_id, -1, 1, '2';
        while (my $seq_id = <MASKED2>) {
            chomp $seq_id;
            $seq_id = substr $seq_id, 1;
            my $mate_seq = <MASKED2>;
            chomp $mate_seq;
            $mate_seq .= <MASKED2>;
            chomp $mate_seq;
            if ($seq_id eq $mate_id) {
                $mate_uniseq = longest_ucstr($mate_seq);
                last;
            }
        }
    } else {
        $mate_id = $id;
        substr $mate_id, -1, 1, '1';
        while (my $seq_id = <MASKED1>) {
            chomp $seq_id;
            $seq_id = substr $seq_id, 1;
            my $mate_seq = <MASKED1>;
            chomp $mate_seq;
            $mate_seq .= <MASKED1>;
            chomp $mate_seq;
            if ($seq_id eq $mate_id) {
                $mate_uniseq = longest_ucstr($mate_seq);
                last;
            }
        }
    }
    
    # evaluate the mappability
    if (length $uniseq < 25) {
        if (length $mate_uniseq < 25) {
            # not mappable
            next;
        } else {
            # mate mappable
            print OUTFILE ">$mate_id\n";
            print OUTFILE "$mate_uniseq\n";
        }        
    } else {
        if (length $mate_uniseq < 25) {
            # mappable, but mate isn't
            print OUTFILE ">$id\n";
            print OUTFILE "$uniseq\n";
        } else {
            # both mappable
            if (length $uniseq < length $mate_uniseq) {
                # mate has longer uniseq
                print OUTFILE ">$mate_id\n";
                print OUTFILE "$mate_uniseq\n";
            } else {
                print OUTFILE ">$id\n";
                print OUTFILE "$uniseq\n";
            }
        }
    }
}
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