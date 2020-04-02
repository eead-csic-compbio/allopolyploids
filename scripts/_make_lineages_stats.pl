#!/usr/bin/perl -w

use strict;
use warnings;
use FindBin '$Bin';
use lib "$Bin";
use polyconfig;

# Takes log made by _check_lineages_polyploids.pl and prints TSV stats
#  
# B Contreras-Moreira, R Sancho, EEAD-CSIC & EPS-UNIZAR 2017-20

my @polyploids = @polyconfig::polyploids;
my @CODES = @polyconfig::CODES; 

my ($taxon,$col,$code,%stats);

die "# usage: $0 <lineage_codes.log>\n" if(!$ARGV[0]);

# read FASTA sequences to get their length
open(LOG,"<",$ARGV[0]);
while(my $line = <LOG>)
{
	#        A       B       C       D       E       F       G       H       I       all
	#Taes    3       1       30      4       1       40      5       5       31      120
	#Ttur    2       1       17      2       1       0       4       7       24      58

	# loop polyploid taxons and match current line	
	foreach $taxon (@polyploids){
		if($line =~ /^$taxon/){
			chomp($line);; 
			my @data = split(/\t/,$line);
			shift(@data);
			foreach $col (0 .. $#data){
				$stats{$taxon}{$CODES[$col]} += $data[$col];
			}
			last;
		} 
	}	
}
close(LOG);

# print stats
foreach $code (@CODES){ print "\t$code" } print "\n";
foreach $taxon (@polyploids){
	print "$taxon";
	foreach $code (@CODES){
		printf("\t%d",$stats{$taxon}{$code} || 0);
	}	print "\n";
}		

