#!/usr/bin/perl -w

# Takes log made by _check_lineages_polyploids.pl and prints TSV stats
#  
# # B Contreras-Moreira, R Sancho, EEAD-CSIC & EPS-UNIZAR 2017

my @polyploids = ('Bhyb','Bboi','Bret','Bmex','Brup','Bpho','B422');
my @CODES = qw( A B C D E F G H I all ); 

my ($taxon,$col,$code,%stats);

die "# usage: $0 <log.lineage_codes>\n" if(!$ARGV[0]);

# read FASTA sequences to get their length
open(LOG,"<",$ARGV[0]);
while(<LOG>)
{
	#	A	B	C	D	E	F	G	H	I	all
	#Bhyb	2	0	0	0	0	0	0	0	0	2
	if(/^(B\S+)/){ 
		$taxon = $1;
		chomp; 
		my @data = split(/\t/,$_);
		shift(@data);
		foreach $col (0 .. $#data){
			$stats{$taxon}{$CODES[$col]} += $data[$col];
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

