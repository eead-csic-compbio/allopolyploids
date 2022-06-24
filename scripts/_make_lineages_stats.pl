#!/usr/bin/perl -w

use strict;
use warnings;
use FindBin '$Bin';
use lib "$Bin";
use polyconfig;

# Takes log made by _check_lineages_polyploids.pl and prints TSV stats,
# which look like this:
#A       B       C       D       E       F       G       H       I       all
#Taes    3       1       30      4       1       40      5       5       31      120
#Ttur    2       1       17      2       1       0       4       7       24      58

# Note: ancestral allele counts with asterisks (ie 1*) can be produced if
# polyconfig::subgenomes are defined. In this case two matrices are produced.

# B Contreras-Moreira, R Sancho, EEAD-CSIC & EPS-UNIZAR 2017-22

my @polyploids = @polyconfig::polyploids;
my @CODES = @polyconfig::CODES; 

my ($taxon,$col,$code,$anc,%stats,%subgenostats);

die "# usage: $0 <lineage_codes.log>\n" if(!$ARGV[0]);

# read FASTA sequences to get their length
open(LOG,"<",$ARGV[0]) || die "# ERROR: cannot read $ARGV[0]\n";
while(my $line = <LOG>) {

  #        A       B       C       D       E       F       H       I       all
  #Bhyb    0       1       0       0       0       0       0       0       1
  #BmeP    1*      0       0       0       0       0       0       0       1
  #BmeU    1       0       0       0       0       0       0       0       1

  # loop polyploid taxons and match current line	
  foreach $taxon (@polyploids){
    if($line =~ /^$taxon/){
      chomp($line);
      my @data = split(/\t/,$line);
      shift(@data);
      foreach $col (0 .. $#data){
        if($data[$col] =~ /(\d+)\*/) {
          $anc = $1;
          $subgenostats{$taxon}{$CODES[$col]} += $anc;
	  $subgenostats{$taxon}{'all'} += $anc;
          $stats{$taxon}{'all'} -= $anc;
        } else {
          $stats{$taxon}{$CODES[$col]} += $data[$col];
        }
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
  } print "\n";
} print "\n";		

# print subgenome stats if required
if(%subgenostats) {
  print "# stats of ancestral subgenome alleles:\n\n";

  foreach $code (@CODES){ print "\t$code" } print "\n";
  foreach $taxon (@polyploids){
    print "$taxon";
    foreach $code (@CODES){
      printf("\t%d",$subgenostats{$taxon}{$code} || 0);
    } print "\n";
  } print "\n";
}
