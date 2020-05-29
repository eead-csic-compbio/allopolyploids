#!/usr/bin/perl

# Re-roots an input Newick tree with an outgroup and 
# prints the resulting sorted tree.
# Outgroup is defined in polyconfig module.
# Based on https://github.com/phac-nml/snvphyl-tools/blob/master/rearrange_snv_matrix.pl
#
# B Contreras-Moreira, R Sancho, EEAD-CSIC & EPS-UNIZAR 2018-20

use strict;
use warnings;
use FindBin '$Bin';
use lib "$Bin";
use polyconfig;

# external dependencies, need to be installed
use Bio::TreeIO;
use Bio::Phylo::IO;
use Bio::Phylo::Forest::Tree;
 
my ($outfound,$outnode,$sorted_newick) = (0);

die "# usage: $0 <tree.newick>\n" if(!$ARGV[0]);

# set outgroup root node
my $outgroup_string = $polyconfig::ROOT;

# read input tree
my $input = new Bio::TreeIO(-file=>$ARGV[0],-format=>'newick');
my $intree = $input->next_tree();

# find outgroup taxon
for my $node ( $intree->get_nodes() ) { 
  if(defined($node->id()) && $node->id() =~ m/$outgroup_string/) {
    $outnode = $node; 
    $outfound = 1;
    last;
  }
}
if($outfound == 0) {
  die "# cannot find outgroup $outgroup_string in input tree $ARGV[0]\n";
}

# root in outgroup and sort in increasing order
my $output = $intree->reroot_at_midpoint($outnode);
if($output == 0) {
	die "# cannot midpoint-root input tree $ARGV[0]\n";
}

# sort nodes in defined order and print
my $unsorted_tree = Bio::Phylo::IO->parse(
  '-string' => $intree->as_text('newick'),
  '-format' => 'newick'
)->first();

$unsorted_tree->ladderize($polyconfig::NODEORDER);

$sorted_newick = $unsorted_tree->to_newick();
$sorted_newick =~ s/'//g;

print $sorted_newick;
