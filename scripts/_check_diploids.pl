#!/usr/bin/perl

# Re-roots an input Newick tree with a user-defined outgroup,
# identifies diploid nodes, prunes the non-diploids and finally 
# outputs the ladderized diploid backbone 
#
# Diploid species must contain unique string hard-coded in @polyconfig::diploids
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

# user-defined regexes to match diploid species
# this names appear appended at end of node names separated by '_'
#my @diploids = ('Tura','Tmon','Asha','Atau','Aspe','Bdis','Hvul','Osat' );
#my %can_be_sisters = ( 'Tura' => 'Tmon', 'Tmon' => 'Tura', 'Asha' => 'Atau', 'Atau' => 'Asha' );

# external dependency, take it from http://cegg.unige.ch/newick_utils 
my $PRUNEXE = '~contrera/soft/newick_utils/src/nw_prune';

######################################################################

my ($outfound,$outnode,$sorted_newick,$taxon,$node,$d,$header) = (0);
my ($non_diploid_labels);

die "# usage: $0 <tree.newick> <outgroup>\n" if(!$ARGV[1]);

my $pruned_tree_file = $ARGV[0].'.pruned';
my $outgroup_string = $ARGV[1];

# read input tree 
my $input = new Bio::TreeIO(
  -file=>$ARGV[0],
  -format=>'newick',
  -internal_node_id => 'bootstrap'
);
my $intree = $input->next_tree();

# check outgroup is there
for $node ( $intree->get_nodes() ) {
  if(defined($node->id()))
  {
    if(!$outfound && $node->id() =~ m/$outgroup_string/) {
      $outnode = $node;
      $outfound = 1;
      last;
    }
  }
}

if($outfound == 0) {
  die "# ERROR: cannot find outgroup $outgroup_string in input tree $ARGV[0]\n";
}

# root in outgroup
$intree->reroot_at_midpoint($outnode);

# copy rooted tree to a new Bio::Phylo::IO object, which supports node laddering
my $sorted_tree = Bio::Phylo::IO->parse(
  '-string' => $intree->as_text('newick'),
  '-format' => 'newick',
  '-internal_node_id' => 'bootstrap'
)->first();
$sorted_tree->ladderize($polyconfig::NODEORDER);

print "# $ARGV[0] $pruned_tree_file\n";

# identify non-diploid labels in sorted tree
for $node ( $sorted_tree->get_nodes() ) {
	next if ($node->id() eq ''); # skip internal nodes

	my $is_diploid = 0; 
	foreach $d ( 0 .. $#polyconfig::diploids) {
		if($node->id() =~ m/($polyconfig::diploids[$d]$)/) {
			$is_diploid = 1;
			last;			
		}
	}

	if($is_diploid == 0) { $non_diploid_labels .= $node->id().' ' }
} 

# prune non diploid leaves
system("$PRUNEXE $ARGV[0] $non_diploid_labels > $pruned_tree_file");

# re-read tree in pruned version
my $pruned = new Bio::TreeIO(
  -file    => $pruned_tree_file,
  -format  => 'newick'
);
my $prunedtree = $pruned->next_tree();

my $sorted_pruned_tree = Bio::Phylo::IO->parse(
  -string => $prunedtree->as_text('newick'),
  -format => 'newick',
)->first();
$sorted_pruned_tree->ladderize($polyconfig::NODEORDER);

# walk the tree from root to tip
my ($total_nodes) = (0);
for $node ( $sorted_pruned_tree->get_nodes() ) {

	if(defined($node->id()) && $node->id() ne '') {
	$total_nodes++;

    # is this a diploid?
    foreach $d ( 0 .. $#polyconfig::diploids) {
      if($node->id() =~ m/($polyconfig::diploids[$d]$)/) {
			$taxon = $1;

			print "# $taxon\t$total_nodes\t".$node->id()."\n";

			last;
      }
    }    
  }
}  

## print simplified Newick topology
my $newick_string = $sorted_pruned_tree->to_newick();

# eliminate distances
$newick_string =~ s/:[\d\.]+//g; 

# simplify node names
foreach $taxon (@diploids){
	$newick_string =~ s/[^\(\),]+?_*$taxon([\(\),])/$taxon$1/g;
}

# terminate by printing out the Newick string 
print "$newick_string\n";
