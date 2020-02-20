#!/usr/bin/perl -w

# Re-roots an input Newick tree with a user-defined outgroup,
# identifies diploid nodes, prunes the non-diploids and finally outputs the ladderized diploid backbone 
#
# Diploid species must contain unique string hard-coded in @diploids
# 
# B Contreras-Moreira, R Sancho, EEAD-CSIC & EPS-UNIZAR 2018-20
# Updated Sep2018 to make sure diploids are not sisters with the exception of %can_be_sisters

use strict;
use Bio::TreeIO;
use Bio::Phylo::IO;
use Bio::Phylo::Forest::Tree;

# user-defined regexes to match diploid species
my @diploids = ('Bdis','Bsta','Barb','Bpin','Bsyl','Hvul','Osat' ); # _Sbic

my %can_be_sisters = ( 'Bpin' => 'Bsyl', 'Bsyl' => 'Bpin' );

# get it from http://cegg.unige.ch/newick_utils
my $PRUNEXE = '/path/to/newick_utils/src/nw_prune';

my $NODEORDER = 0; # 1:increasing, 0:decreasing
my ($outfound,$outnode,$sorted_newick,$taxon,$node,$d,$header) = (0);
my ($previous_node,$previous_taxon,$is_diploid,$non_diploid_labels);
my ($bases,%length,%longest,%longest_header);

die "# usage: $0 <tree.newick> <outgroup>\n" if(!$ARGV[1]);

# name of output pruned file
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
  die "# cannot find outgroup $outgroup_string in input tree $ARGV[0]\n";
}

# root in outgroup
$intree->reroot_at_midpoint($outnode);

# copy rooted tree to a new Bio::Phylo::IO object, which supports node laddering
my $sorted_tree = Bio::Phylo::IO->parse(
  '-string' => $intree->as_text('newick'),
  '-format' => 'newick',
  '-internal_node_id' => 'bootstrap'
)->first();
$sorted_tree->ladderize($NODEORDER);

print "# $ARGV[0] $pruned_tree_file\n";

# identify non-diploid labels in sorted tree
for $node ( $sorted_tree->get_nodes() ) {
	next if ($node->id() eq ''); # skip internal nodes

	$is_diploid = 0;
	foreach $d ( 0 .. $#diploids) {
		if($node->id() =~ m/_($diploids[$d])/) {
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
$sorted_pruned_tree->ladderize($NODEORDER);

# walk the tree from root to tip
my ($total_nodes,$summary,%seen) = (0,'');
$previous_node = undef;
$previous_taxon = '';
LOOP: for $node ( $sorted_pruned_tree->get_nodes() ) {

# debugging: make sure we understand internal ids and node ancestors
#$taxon = 'NA';
#foreach $d ( 0 .. $#diploids) {
#		if($node->id() ne '' && $node->id() =~ m/_($diploids[$d])/) {
#         $taxon = $1;
#			last;
#		}
#	}
#	my @these_ancestors = get_all_ancestors( $node );	
#	printf("# node %s = %s %s ancestry: %s\n",$node->internal_id(),$node->id() || 'internal',$taxon,join(',',@these_ancestors)); # debug	

	if(defined($node->id()) && $node->id() ne '') {
	$total_nodes++;

    # is this a diploid?
    foreach $d ( 0 .. $#diploids) {
      if($node->id() =~ m/_($diploids[$d])/) {
			$taxon = $1;

			print "# $taxon\t$total_nodes\t".$node->id()."\n";

			# check whether this diploid is sister to previous diploid in ladderized tree
			if(defined($previous_node) && 
				(!$can_be_sisters{$taxon} || $can_be_sisters{$taxon} ne $previous_taxon)) {

				#print "check sisters $previous_taxon <-> $taxon\n"; # debug

				my @ancestors_previous_taxon = get_all_ancestors( $previous_node );
				my @ancestors_taxon = get_all_ancestors( $node );

				#my $MRCA = get_MRCA( \@ancestors_previous_taxon, \@ancestors_taxon );
				#print $previous_node->ancestor()->internal_id()." ".$node->ancestor()->internal_id()." $MRCA\n"; # debug

				if($previous_node->ancestor()->internal_id() eq $node->ancestor()->internal_id() ){
					# nodes share immediate ancestor
					$summary .= "$taxon*,";				
				}
				else {
					$summary .= "$taxon,";
				}
			}
			else {
				$summary .= "$taxon,";
			}

			$previous_node = $node;
			$previous_taxon = $taxon;
			last;
      }
    }    
  }
}  

print "$summary\n";

# get recursively all ancestors of a node 
# returns a list of internal ids sorted from immediate to toplevel/root ancestor
sub get_all_ancestors
{
   my ($node,$verbose) = (@_);
   my @ancestors;

	

   while(defined($node->ancestor()))
   {
      $node = $node->ancestor();
		push(@ancestors,$node->internal_id());
   }

	return @ancestors;
}

# receives two sorted lists of ancestors 
# returns the internal_id of MRCA (most_recent_common_ancestor)
sub get_MRCA
{
   my ($ref_ancestorsA, $ref_ancestorsB) = @_;
   my ($nodeA,$nodeB,%hashB);
   my $MRCA = undef;

	# hash elements in list B
	foreach $nodeB (@$ref_ancestorsB)
   {
      $hashB{$nodeB} = 1;
   }

	# iterate along A
	foreach $nodeA (@$ref_ancestorsA)
   {
      if($hashB{$nodeA})
      {
         $MRCA = $nodeA;
         last;
      }
   }

   return $MRCA;
}
