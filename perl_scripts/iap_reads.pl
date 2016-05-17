
#!/usr/bin/perl -w
use strict;

# Screen the RepeatMasker output for reads containing IAP
# Author: Ying Zhang; yzhang@bccrc.ca
# Date: Nov 27, 2014

my $dataset = "c3h.input.2";
my $dir = "$ENV{HOME}/Work/BCCRC/Projects/c3h_chipseq/data";
my $output_file = "$dir/rm_out/$dataset/$dataset.iap";

open (OUTFILE, ">$output_file") or die "cannot open OUTPUT file: $!";

my %iap;
my $total_files = 600;
for (my $i=1; $i<=$total_files; $i++) {
    my $rm_file = "$dir/rm_out/$dataset/$dataset.part$i.fa.out";
    open (RMFILE, "$rm_file") or die "cannot open OUTPUT file: $!";
    <RMFILE>;
    <RMFILE>;
    <RMFILE>;
    
    while (my $input_line = <RMFILE>) {
        if ($input_line =~ /IAP/) {
            my @fields = split /\s+/, $input_line;
            my $id = $fields[5];
            $id = substr $id, 0, -2;
            @fields = split /:/, $id;
            $iap{$fields[1]}{$fields[2]}{$fields[3]}{$fields[4]} = $input_line;
        }
    }
    close RMFILE;
}

foreach my $iap1 (sort {$a <=> $b} keys %iap) {
    foreach my $iap2 (sort {$a <=> $b} keys $iap{$iap1}) {
        foreach my $iap3 (sort {$a <=> $b} keys $iap{$iap1}{$iap2}) {
            foreach my $iap4 (sort {$a <=> $b} keys $iap{$iap1}{$iap2}{$iap3}) {
                print OUTFILE $iap{$iap1}{$iap2}{$iap3}{$iap4};
            }
        }
    }
}

print "Done!\n";
