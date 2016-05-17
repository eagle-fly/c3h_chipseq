
#!/usr/bin/perl -w
use strict;

# Remove regions that are overlapping in bedGraph file
# Author: Ying Zhang; yzhang@bccrc.ca
# Date: Feb 12, 2015

my $dataset = "129.k4m3";
my $dir = "$ENV{HOME}/Work/BCCRC/Projects/c3h_chipseq/results";
my $input_file ="$dir/$dataset"."_treat_pileup_mm9sorted.bdg";
open (INFILE, "$input_file") or die "cannot open INPUT file: $!";
my $output_file = "$dir/$dataset"."_treat_pileup_mm9clean.bdg";
open (OUTFILE, ">$output_file") or die "cannot open OUTPUT file: $!";

my $input_line = <INFILE>;
chomp $input_line;
my ($chr0, $start0, $end0, $value0) = split /\t/, $input_line;

my $overlapped=0;
while ($input_line = <INFILE>) {
    chomp $input_line;
    my ($chr, $start, $end, $value) = split /\t/, $input_line;
    
    # check if overlap
    if ($chr eq $chr0 and $start < $end0) {
        $overlapped = 1;
    } else {
        print OUTFILE "$chr0\t$start0\t$end0\t$value0\n";
        $overlapped = 0;
    }
    
    $chr0 = $chr;
    $start0 = $start;
    $end0 = $end;
    $value0 = $value;
}
# check if the last one is not overlapped
if ($overlapped == 0) {
    print OUTFILE "$chr0\t$start0\t$end0\t$value0\n";
}

print "Done!\n";
