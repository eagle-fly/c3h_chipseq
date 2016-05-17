
#!/usr/bin/perl -w
use strict;

# Merge masked reads into a single file
# Author: Ying Zhang; yzhang@bccrc.ca
# Date: Nov 27, 2014

my $dataset = "c3h.input.2";
my $output_file = "../data/rm_out/$dataset/$dataset.fa.masked";

open (OUTFILE, ">$output_file") or die "cannot open OUTPUT file: $!";

my %masked;
my $total_files = 600;
for (my $i=1; $i<=$total_files; $i++) {
    my $masked_file = "../data/rm_out/$dataset/$dataset.part$i.fa.masked";
    open (MASKED, "$masked_file") or die "cannot open masked file: $!";
    while (my $input_line = <MASKED>) {
        my $id = $input_line;
        chomp $id;
        $id = substr $id, 0, -2;
        my @fields = split /:/, $id;
        my $seq = <MASKED>;
        $seq .= <MASKED>;
        $masked{$fields[1]}{$fields[2]}{$fields[3]}{$fields[4]} = $input_line.$seq;
    }
    close MASKED;
}

foreach my $masked1 (sort {$a <=> $b} keys %masked) {
    foreach my $masked2 (sort {$a <=> $b} keys $masked{$masked1}) {
        foreach my $masked3 (sort {$a <=> $b} keys $masked{$masked1}{$masked2}) {
            foreach my $masked4 (sort {$a <=> $b} keys $masked{$masked1}{$masked2}{$masked3}) {
                print OUTFILE $masked{$masked1}{$masked2}{$masked3}{$masked4};
            }
        }
    }
}

print "Done!\n";
