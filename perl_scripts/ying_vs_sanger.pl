
#!/usr/bin/perl -w
use strict;

# Compare C3H-specific IAPs between Ying's data and Sanger's data
# Author: Ying Zhang; yzhang@bccrc.ca
# Date: Jan 15, 2015


use DBI;

   
#---------------------------------------------------------------------------------------
# Perform the connection using the mysql driver,
# then construct the string to pass to DBI
#---------------------------------------------------------------------------------------

my $database = "polyte";
my $dbh = DBI->connect( "dbi:mysqlPP:$database", "root", "bccrc", 
	{
		PrintError => 0,
            	RaiseError => 0 
    	} 
) or die "Can't connect to the database: $DBI::errstr\n";

#---------------------------------------------------------------------------------------

my $ying_data = "../results/uniseq/c3h.iap.bed";
open (YINGDATA, "$ying_data") or die "cannot open Ying's Data: $!";
my $output_file = "../results/uniseq/ying_vs_sanger.txt";
open (OUTPUT, ">$output_file") or die "cannot open output file: $!";
print OUTPUT "chrom\tstart\tend\tying\tsanger\n";

<YINGDATA>;
while (<YINGDATA>) {
    chomp;
    my ($chr, $start, $end, @fields) = split /\t/;
    my $last_column = $fields[-1];
    my $status;
    my $offset;
    if ($last_column != 0) {
	$status = "c3h_only";
	$offset = 200;
    } else {
	$status = "common";
	$offset = 10000;
    }
        
    # check if the insertion exists in the Sanger dataset
    my $sql = "SELECT chrom, start, end FROM C3H_polyIAP_Sanger
                WHERE chrom = \"$chr\"
                                        AND ((start >= $start - $offset AND start <= $start + $offset)
                                        OR (end >= $end - $offset AND end <= $end + $offset))";
    my $sth = $dbh->prepare($sql) or die "Can't prepare SQL statement: $DBI::errstr\n";
    $sth->execute or die "Can't execute SQL statement: $DBI::errstr\n";
    my @rows = $sth->fetchrow_array;
    my $sanger = "";
    if ($rows[0] ne "") {
	$sanger = "c3h_only";
    }
    
    print OUTPUT "$chr\t$start\t$end\t$status\t$sanger\n";
}

print "Done!\n";
