
#!/usr/bin/perl -w
use strict;

# convert the track_data file from mm9 to mm10 genome
# Author: Ying Zhang; yzhang@bccrc.ca
# Date: Oct 23, 2014

my $old_track = '../results/tt2.k4m3_peaks.narrowPeak';
my $old_loci = '../results/track_loci_mm10.txt';
my $new_loci = '../results/track_loci_mm9.txt';
my $new_track = '../results/tt2.k4m3_peaks_mm9.narrowPeak';

open (OLDTRACK, "$old_track") or die "cannot open track data file: $!";
open (OLDLOCI, ">$old_loci") or die "cannot open track data file: $!";

my %data;
while (<OLDTRACK>) {
    if (!/^chr/) {
        next;
    }
    my ($chr, $start, $end, $name, @other_fields) = split /\t/, $_;
    print OLDLOCI "$chr\t$start\t$end\t$name\n";
    # store data into a hash
    $data{$name} = join "\t", @other_fields;
}

close OLDTRACK;
close OLDLOCI;

# convert old loci from mm10 to mm9 by using the liftOver tool
my $liftover = "./liftOver ../results/track_loci_mm10.txt mm10ToMm9.over.chain ../results/track_loci_mm9.txt ../results/mm9.unmaped.txt";
system ($liftover) and die "cannot liftOver\n";

open (OLDTRACK, "$old_track") or die "cannot open track data file: $!";
open (NEWLOCI, "$new_loci") or die "cannot open track data file: $!";
open (NEWTRACK, ">$new_track") or die "cannot open track description file: $!";

# print track configuration lines
my $line = <OLDTRACK>;
do {
    $line =~ s/mm10/mm9/;
    $line =~ s/k4m3_peaks.narrowPeak/k4m3_peaks_mm9.narrowPeak/;
    print NEWTRACK $line;
    $line = <OLDTRACK>;
} until ($line =~ /^chr/);

# print data lines with new loci info
while (<NEWLOCI>) {
    chomp;
    my ($chr, $start, $end, $name) = split /\t/, $_;
    print NEWTRACK "$chr\t$start\t$end\t$name\t$data{$name}";
}

close OLDTRACK;
close NEWLOCI;
close NEWTRACK;

print "\nDone!\n";