
#!/usr/bin/perl -w
use strict;

# Split the original FASTA files into smaller files before RepeatMasker
# and create the corresponding job batch files
# Author: Ying Zhang; yzhang@bccrc.ca
# Date: Nov 24, 2014

my $dataset = "c3h.k4m3.2";
my $input_file = "../data/mapping/$dataset.fa";
open (INFILE, "$input_file") or die "cannot open input file: $!";
my $qsub_file = "../data/rm_out/qsub_$dataset.sh";
open (QSUBFILE, ">$qsub_file") or die "cannot open OUTPUT file: $!";

my $total_lines = 90711114;
my $total_files = 900;
my $line_max = int($total_lines / $total_files);
for (my $i=1; $i<$total_files; $i++) {
    my $output_file = "../data/split/$dataset.part$i.fa";
    open (OUTFILE, ">$output_file") or die "cannot open OUTPUT file: $!";
    for (my $line=1; $line<=$line_max; $line+=2) {
        my $read_id = <INFILE>;
        print OUTFILE $read_id;
        my $read_seq = <INFILE>;
        print OUTFILE $read_seq;
    }
    close OUTFILE;
    my $batch_file = "../data/rm_out/rm.$dataset.part$i.sh";
    open (BATCHFILE, ">$batch_file") or die "cannot open batch file: $!";
    print BATCHFILE "#\!/bin/sh\n";
    print BATCHFILE "#\$ -S /bin/sh\n";
    print BATCHFILE "#\$ -m e\n";
    print BATCHFILE "#\$ -M yzhang\@bccrc.ca\n";
    print BATCHFILE "#\$ -N RM$i\n";
    print BATCHFILE "#\$ -V\n";
    print BATCHFILE "RepeatMasker -xsmall -species mouse -no_is -dir $dataset ../split/$dataset.part$i.fa\n";
    close BATCHFILE;
    
    print QSUBFILE "qsub -l mem_free=8G -l mem_token=8G rm.$dataset.part$i.sh\n";
}
# last file
my $output_file = "../data/split/$dataset.part$total_files.fa";
open (OUTFILE, ">$output_file") or die "cannot open OUTPUT file: $!";
while (<INFILE>) {
    print OUTFILE;
}
my $batch_file = "../data/rm_out/rm.$dataset.part$total_files.sh";
open (BATCHFILE, ">$batch_file") or die "cannot open batch file: $!";
print BATCHFILE "#\!/bin/sh\n";
print BATCHFILE "#\$ -S /bin/sh\n";
print BATCHFILE "#\$ -m e\n";
print BATCHFILE "#\$ -N RM$total_files\n";
print BATCHFILE "#\$ -V\n";
print BATCHFILE "RepeatMasker -xsmall -species mouse -no_is -dir $dataset ../split/$dataset.part$total_files.fa\n";
close BATCHFILE;
print QSUBFILE "qsub -l mem_free=8G -l mem_token=8G rm.$dataset.part$total_files.sh\n";

print "Done!\n";
