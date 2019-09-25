#!/usr/bin/perl -w

# Re-roots an input Newick tree with a user-defined outgroup and 
# prints the resulting tree in ascending or descending node order.
# Based on https://github.com/phac-nml/snvphyl-tools/blob/master/rearrange_snv_matrix.pl
#
# B Contreras-Moreira, R Sancho, EEAD-CSIC & EPS-UNIZAR 2018


use strict;
use Bio::TreeIO;
use Bio::Phylo::IO;
use Bio::Phylo::Forest::Tree;
 
my $OUTGROUPSTRING = 'Osat'; # change as needed _Sbic
my $NODEORDER = 0; # 1:increasing, 0:decreasing
my ($outfound,$outnode,$sorted_newick) = (0);

die "# usage: $0 <tree.newick>\n" if(!$ARGV[0]);

# read input tree
my $input = new Bio::TreeIO(-file=>$ARGV[0],-format=>'newick');
my $intree = $input->next_tree();

# find outgroup taxon
for my $node ( $intree->get_nodes() ) { 
  if(defined($node->id()) && $node->id() =~ m/$OUTGROUPSTRING/) {
    $outnode = $node; 
    $outfound = 1;
    last;
  }
}
if($outfound == 0) {
  die "# cannot find outgroup $OUTGROUPSTRING in input tree $ARGV[0]\n";
}

# root in outgroup and sort in increasing order
$intree->reroot_at_midpoint($outnode);

# sort nodes in defined order and print
my $unsorted_tree = Bio::Phylo::IO->parse(
  '-string' => $intree->as_text('newick'),
  '-format' => 'newick'
)->first();

$unsorted_tree->ladderize($NODEORDER);

$sorted_newick = $unsorted_tree->to_newick();
$sorted_newick =~ s/'//g;

print $sorted_newick;
