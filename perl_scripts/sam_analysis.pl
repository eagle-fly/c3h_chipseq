
#!/usr/bin/perl -w
use strict;

# parse the SAM file for mapping positions
# and consolidate overlapped positions
# Author: Ying Zhang; yzhang@bccrc.ca
# Date: Dec 1, 2014


use DBI;

   
#---------------------------------------------------------------------------------------
# Perform the connection using the mysql driver,
# then construct the string to pass to DBI
#---------------------------------------------------------------------------------------

my $database = "rmsk";
my $dbh = DBI->connect( "dbi:mysqlPP:$database", "root", "bccrc", 
	{
		PrintError => 0,
            	RaiseError => 0 
    	} 
) or die "Can't connect to the database: $DBI::errstr\n";

#---------------------------------------------------------------------------------------

my $dataset = "tt2";
my $sam_file = "../results/uniseq/$dataset.sam";
open (SAMFILE, "$sam_file") or die "cannot open sam file: $!";
my $bed_file = "../results/uniseq/$dataset.iap2.bed";
open (BEDFILE, ">$bed_file") or die "cannot open BED file: $!";
print BEDFILE "track type=bed itemRgb=On visibility=full db=mm10 name=\"IAPs in C3H\" description=\"IAPs detected using the pooled C3H ChIP-seq reads\"\n";

my %iaps;
while (<SAMFILE>) {
    # skip the header lines
    if (/^@/) {
        next;
    }
    my @fields = split /\t/;
    my $read_id = $fields[0];
    my $chr = $fields[2];
    my $pos = $fields[3];
    my $cigar = $fields[5];

    if ($chr eq "*") { # unmapped
        next;
    }
    
    $iaps{$chr}{$pos} = $read_id;    
}

my $last_chr = "chr1";
my $last_pos = 0;
my $read_cutoff = 2;
my $read_support = 1;
foreach my $chr (sort keys %iaps) {
    foreach my $pos (sort keys $iaps{$chr}) {
        if ($chr eq $last_chr and abs($pos - $last_pos) < 10000) { # overlapping position; skip
            delete $iaps{$last_chr}{$last_pos};
            $read_support++;
        } else {
            if ($read_support < $read_cutoff) {
                delete $iaps{$last_chr}{$last_pos};
            }
            $read_support = 1;
        }
        $last_chr = $chr;
        $last_pos = $pos;
    }
}

foreach my $chr (sort keys %iaps) {
    foreach my $pos (sort keys $iaps{$chr}) {
        my $start = $pos - 200;
        my $end = $pos + 200;
        my $name = $iaps{$chr}{$pos};
        my $overlap = iap_overlap($chr, $start, $end);
        print BEDFILE "$chr\t$start\t$end\t$name\t999\t+\t$start\t$end\t$overlap\n";
    }
}
print "Done!\n";

sub iap_overlap {
    my $chr = shift;
    my $iStart = shift;
    my $iEnd = shift;
    
    # check if the IAP insertion in C3H strain is overlapping with those in B6
    my $sql = "SELECT genoName, genoStart, genoEnd FROM mm10_rmsk_iap
                WHERE genoName = \"$chr\"
                                        AND ((genoStart <= $iStart AND genoEnd >= $iStart)
                                        OR (genoStart <= $iEnd AND genoEnd >= $iEnd)
                                        OR (genoStart >= $iStart AND genoEnd <= $iEnd))";
    my $sth = $dbh->prepare($sql) or die "Can't prepare SQL statement: $DBI::errstr\n";
    $sth->execute or die "Can't execute SQL statement: $DBI::errstr\n";
    my @rows = $sth->fetchrow_array;
    my $start = $rows[1];
    my $end = $rows[2];
    if ($start ne "") {
        return 0
    } else {
        return "255,0,0"
    }
} 