#!/usr/bin/perl -w

use strict;
use IO::String;
use Bio::TreeIO;
use Bio::Phylo::IO;
use Bio::Phylo::Forest::Tree;

# Label allopolyploid sequences with @CODES according to their position
# with respect to the diploid backbone. When several species of the same
# allopolyploid species are found, all of them are labelled.
#  
# INPUT: Reads multiple sequence alignment (MSA) corresponding to input tree file (.root.ph) 
# OUTPUT: produces two labelled MSA files with sorted diploid and polyploid (A,B,C,etc) species,
# a full MSA with some only-gap rows, and a reduced MSA excluding those.
# Ad-hoc rules encoded herein (see below) define pairs of diploid species where one on the can be 
# missing in the alignment. In this case they are Osat/Hvul and Bpin/Bsyl
#
# B Contreras-Moreira, R Sancho, EEAD-CSIC & EPS-UNIZAR september 2018

my $DEBUG = 1;

my @diploids = qw( Osat Hvul Bsta Bdis Barb Bpin Bsyl ); # expects one seq per species and tree
my @polyploids = ('Bhyb','Bboi','Bret','Bmex','Brup','Bpho','B422');
my $NODEORDER = 0; # 1:increasing, 0:decreasing

my @CODES = qw( A B C D E F G H I all ); # see hard-coded rules below

my ($outfound,$outnode,$sorted_newick,$taxon,$node,$node_id,$child,$anc) = (0);
my ($min_anc_dist,$min_desc_dist,$dist,$anc_dip_taxon,$desc_dip_taxon) = (0,0);
my ($coord,$coord5,$coord3,$midcoord);
my $MSAwidth = 0;
my ($d,$header,$node1,$node2,$taxon1,$taxon2,$bases,$lineage_code,$full_taxon);
my ($poly_ancestor_diploid_MRCA,$desc_is_sister);
my (%length,%FASTA,%order,%dip_taxon,%dip_MRCA,%MRCA_computed,%dip_ancestor,%header2label);
my (%dip_all_ancestors,%dip_all_ancestors_string,$MRCA,$diploid_MRCA,$diploid_MRCA_node);

die "# usage: $0 <treefile.root.ph>\n" if(!$ARGV[0]);

my $tree_file = $ARGV[0];
my $MSA_file = $tree_file;
$MSA_file =~ s/\.treefile\.root\.ph//;
my $labelledMSA_file = $MSA_file;
$labelledMSA_file =~ s/\.fna/.label.fna/;
my $reducedMSA_file = $MSA_file;
$reducedMSA_file =~ s/\.fna/.label.reduced.fna/;
my $labelled_tree_file = $tree_file;
$labelled_tree_file =~ s/\.ph/.label.ph/;

# read FASTA sequences and record length
open(MSA,"<",$MSA_file);
while(<MSA>)
{
	if(/^>(\S+)/){ 
		$header = $1; #print "$header\n";
		$MSAwidth = 0;
	}	
	else{ 
		chomp;
		$FASTA{$header} .= $_;
		$bases = ($_ =~ tr/[ACGTN]//);
		$length{$header} += $bases; 
		$MSAwidth += length($_);
	}	
}
close(MSA);

# ladderize/sort input tree 
my $unsorted_input_tree = Bio::Phylo::IO->parse(
  '-file' => $ARGV[0],
  '-format' => 'newick'
)->first();

$unsorted_input_tree->ladderize($NODEORDER);
$sorted_newick = $unsorted_input_tree->to_newick();
$sorted_newick =~ s/'//g;
my $iostring = IO::String->new($sorted_newick);


# read pre-sorted input tree
my $input = new Bio::TreeIO(
  -fh=>$iostring,
  -format=>'newick',
);
my $intree = $input->next_tree();


# find reference diploid sequences
my @refDiploid;
my $total_nodes = 0;
for $node ( $intree->get_nodes() ) {
	if(defined($node->id())) {
		$total_nodes++;
		foreach $d ( 0 .. $#diploids) {
			if($node->id() =~ m/_($diploids[$d])$/) { 
				push(@refDiploid,$node); 
				$dip_taxon{$node} = $1;
				$order{$node->id()} = $total_nodes; 
				$dip_ancestor{$node} = $node->ancestor(); # immediate ancestor
				print "$total_nodes ".$node->id()." ".$node->ancestor()."\n" if($DEBUG);

				my @all_ancestors = get_all_ancestors($node);
				$dip_all_ancestors_string{$node} = join(',',@all_ancestors);
				$dip_all_ancestors{$node} = \@all_ancestors;  
				print "$dip_all_ancestors_string{$node}\n" if($DEBUG);
			}
		}
	}	
}


# identify most recent common ancestors among diploid nodes
# diploids nodes are compared pairwise in the order they appear from root to tips

foreach $node (0 .. $#refDiploid-1){

	$node1 = $refDiploid[$node];
	$node2 = $refDiploid[$node+1]; 

   if(!defined($MRCA_computed{$node1})){ # sister of previous node

		# assumes node 1 is ancestor to node 2, note that $MRCA is an internal node id, which is an arbitrary integer
		$MRCA = get_MRCA( $dip_all_ancestors{$node1} , $dip_all_ancestors{$node2} ); 
		$dip_MRCA{ $MRCA } = $dip_taxon{$node1};
		$MRCA_computed{$node1} = 1;
		print "MRCA $dip_MRCA{ $MRCA } $MRCA\n" if($DEBUG); 

		# but we should check whether they are sisters
		if($dip_all_ancestors_string{$node1} eq $dip_all_ancestors_string{$node2}){
			print "MRCA $dip_taxon{ $node2 } $MRCA (sister)\n" if($DEBUG);
			$MRCA_computed{$node2} = 1;
		}	
   }

	# add also ancestor of last diploid in order
	if($node == $#refDiploid-1){
		$MRCA = $node2->ancestor()->internal_id(); # fake MRCA for last diploid
		$dip_MRCA{ $MRCA } = $dip_taxon{$node2};
      $MRCA_computed{$node2} = 1;
		print "MRCA $dip_MRCA{ $MRCA } $MRCA\n" if($DEBUG); 
	}
}


# create labelled MSA files
open(LABELMSA,">",$labelledMSA_file) || die "# ERROR: cannot create $labelledMSA_file\n";

open(LABELMSAREDUCED,">",$reducedMSA_file) || die "# ERROR: cannot create $reducedMSA_file\n";

# write output file names
print "\n# $labelledMSA_file $reducedMSA_file $MSAwidth $labelled_tree_file\n\n";

# print diploid sequences
foreach $taxon (@diploids){
	my $matched = 0;
	foreach $node (@refDiploid){
		$header = $node->id();
		if($header =~ $taxon){			 
			
			printf(LABELMSA ">%s %s\n%s\n",$taxon,$header,$FASTA{$header}); 

			printf(LABELMSAREDUCED ">%s %s\n%s\n",$taxon,$header,$FASTA{$header});
			$matched = 1;
         last;
		}            
	}

	# add non-matched diploids for consistency to allow posterior concat
	if(!$matched)
   {
      printf(LABELMSA ">%s\n%s\n",$taxon,'-' x $MSAwidth);
   }
}

# print header of table of labels (see @CODES)
foreach $lineage_code (@CODES){ printf("\t%s",$lineage_code) }  print "\n";

# find polyploid nodes, label nr and add them to labelled MSA
foreach $taxon (@polyploids)
{
	my (%taxon_ancestors,@nr_taxon_nodes,%taxon_order,%taxon_stats,%taxon_seqs);

	print "> $taxon\n" if($DEBUG);

	$total_nodes = 0;
	for $node ( $intree->get_nodes() ) {
		if(defined($node->id())) {
			$total_nodes++;
			if($node->id() =~ m/$taxon$/){

				# found a polyploid sequence			
				$taxon_order{$node->id()} = $total_nodes;
				$anc = $node->ancestor();
				push(@{$taxon_ancestors{$anc}},$node);
			}
		}  
	}

	# simplify same-taxon polytomies by taking longest sequence
	foreach $anc (keys(%taxon_ancestors)){
		my $total_polytomic_nodes = 0;
		foreach $node (sort {$length{$b->id()}<=>$length{$a->id()}} @{$taxon_ancestors{$anc}}){
			$total_polytomic_nodes++;
			if($total_polytomic_nodes == 1){
				push(@nr_taxon_nodes,$node); #printf("%s %d\n",$node->id(),$length{$node->id()});
			}
		}
	}

	# check closest diploid ancestor and descendant of this polyploid node
	foreach $node1 (@nr_taxon_nodes){

		($anc_dip_taxon,$desc_dip_taxon,$desc_is_sister,$lineage_code) = ('-','-',0,'-');

		# get all ancestors of this polyploid node, from immediate to toplevel/root
		my @node1_all_ancestors = get_all_ancestors($node1); 
		my $node1_all_ancestors_string = join(',',@node1_all_ancestors);

		print $node1->id()." $node1_all_ancestors_string\n" if($DEBUG); 

		# find closest ancestor which is MRCA of two diploids (or fake MRCA of last diploid)
		$poly_ancestor_diploid_MRCA = -999;
		foreach $node_id (@node1_all_ancestors){
			if(defined($dip_MRCA{ $node_id })){
				$poly_ancestor_diploid_MRCA = $node_id;
				$anc_dip_taxon = $dip_MRCA{ $node_id }; 
				print "anc_dip_taxon $anc_dip_taxon $node_id\n" if($DEBUG); 
				last;
			}
		}
	
		# find the first descendant which is either i) a diploid sister or ii) a MRCA between diploids
		my %total_shared_desc; 
		foreach $node2 (@refDiploid){
		
			# if all ancestors are shared, then polyploid (node1) and diploid (node2) are sisters			
			if($node1_all_ancestors_string eq $dip_all_ancestors_string{$node2}){
				$desc_is_sister = 1;
				$desc_dip_taxon = $dip_taxon{ $node2 };
				$total_shared_desc{$desc_dip_taxon} = scalar(@{$dip_all_ancestors{$node2}});
				last;
				#print "des_tip_taxon $desc_dip_taxon\n" if($DEBUG); 
			}
			else{ # check shared ancestors between to consecutive diploid MRCA nodes

				my @shared_nodes = get_shared_ancestors_next_diploid_MRCA( $dip_all_ancestors{$node2}, 
					\@node1_all_ancestors, $poly_ancestor_diploid_MRCA, \%dip_MRCA );

				print "<$dip_taxon{ $node2 }> ".join(',',@shared_nodes)." ($poly_ancestor_diploid_MRCA) ". 
					join(',',@{$dip_all_ancestors{$node2}})." ".join(',',@node1_all_ancestors)."\n" if($DEBUG);

				if(scalar(@shared_nodes) > 1){
								
					if($dip_MRCA{ $shared_nodes[0] }){ 
						$desc_dip_taxon = $dip_MRCA{ $shared_nodes[0] };
						$total_shared_desc{$desc_dip_taxon} = scalar(@shared_nodes)-1;
						#print "des_dip_taxon $desc_dip_taxon $shared_nodes[0]\n" if($DEBUG);
					}
					else{
						$desc_is_sister = 1;
		            $desc_dip_taxon = $dip_taxon{ $node2 };
						$total_shared_desc{$desc_dip_taxon} = 1;
				      #print "des_tip_taxon $desc_dip_taxon\n" if($DEBUG);
					}
				}
			}

			#last if($desc_dip_taxon ne '-');
		}

		# check diploid with most shared ancestors, which will be descendant
		my $max_shared = 0;
		$desc_dip_taxon = '';
		foreach $taxon (keys(%total_shared_desc)){
			if($total_shared_desc{$taxon} > $max_shared){
				$desc_dip_taxon = $taxon;
				$max_shared = $total_shared_desc{$taxon};
			}
			#print "> $taxon $total_shared_desc{$taxon}\n";
		}
		print "des_tip_taxon $desc_dip_taxon\n" if($DEBUG);
		#exit;

      # rules based on node order (NEW)
      if($anc_dip_taxon eq 'Hvul'){
         if($desc_dip_taxon eq 'Bsta') {$lineage_code = 'A'}
		}
		elsif($anc_dip_taxon eq 'Bsta'){
			if($desc_dip_taxon eq 'Bdis'){ $lineage_code = 'C' }
			elsif($desc_dip_taxon eq 'Bsta' && $desc_is_sister == 1){ $lineage_code = 'B' }
			#elsif($desc_dip_taxon eq 'Barb'){ $lineage_code = 'B' } # if Bdis and Barb are flipped
			else{ $lineage_code = 'C' }
		}
		elsif($anc_dip_taxon eq 'Bdis'){
			if($desc_dip_taxon eq 'Barb'){ $lineage_code = 'E' }
			elsif($desc_dip_taxon eq 'Bdis' && $desc_is_sister == 1){ $lineage_code = 'D' }
			elsif($desc_dip_taxon eq 'Bpin' || $desc_dip_taxon eq 'Bsyl'){ $lineage_code = 'G' } # if Bpin/Bsyl and Barb are flipped
			else{ $lineage_code = 'D' }
		}
		elsif($anc_dip_taxon eq 'Barb'){
         if($desc_dip_taxon eq 'Bpin' || $desc_dip_taxon eq 'Bsyl'){ $lineage_code = 'G' }
			elsif($desc_dip_taxon eq 'Barb' && $desc_is_sister == 1){ $lineage_code = 'F' }
			else{ $lineage_code = 'F' } # if Bpin/Bsyl and Barb are flipped and Barb are the last diploid species
		}
      elsif($anc_dip_taxon eq 'Bsyl'){
         if($desc_dip_taxon eq 'Bpin'){ $lineage_code = 'H' }
         elsif($desc_dip_taxon eq 'Bsyl' && $desc_is_sister == 1){ $lineage_code = 'H' }
			else{ $lineage_code = 'H' }
      }
		elsif($anc_dip_taxon eq 'Bpin'){
         if($desc_dip_taxon eq 'Bsyl'){ $lineage_code = 'I' }
			elsif($desc_dip_taxon eq 'Bpin' && $desc_is_sister == 1){ $lineage_code = 'I' }
			else{ $lineage_code = 'I' }
		}
		#RUBEN (esto no esta funcionando)
		elsif(($anc_dip_taxon eq 'Bsyl' || $anc_dip_taxon eq 'Bpin') && $desc_dip_taxon eq 'Barb'){
			if($desc_dip_taxon eq 'Bpin' && $desc_is_sister == 1){ $lineage_code = 'I' }	
			elsif($desc_dip_taxon eq 'Bsyl' && $desc_is_sister == 1){ $lineage_code = 'H' }
			else{$lineage_code = 'F'}
		}








		# rules based on node order (ORIGINAL RULES)
#		if(($anc_dip_taxon eq 'Osat' || $anc_dip_taxon eq 'Hvul') && $desc_dip_taxon eq 'Bsta'){
#			$lineage_code = 'A'; 
#		}
#		elsif($anc_dip_taxon eq 'Bsta'){
#			if($desc_dip_taxon eq 'Bdis'){ $lineage_code = 'C' }
#			elsif($desc_dip_taxon eq 'Bsta' && $desc_is_sister == 1){ $lineage_code = 'B' }
#		}
#		elsif($anc_dip_taxon eq 'Bdis'){
#			if($desc_dip_taxon eq 'Barb'){ $lineage_code = 'E' }
#			elsif($desc_dip_taxon eq 'Bdis' && $desc_is_sister == 1){ $lineage_code = 'D' }
#      }
#		elsif($anc_dip_taxon eq 'Barb'){
#			if($desc_dip_taxon eq 'Barb' && $desc_is_sister == 1){ $lineage_code = 'F' }
#			elsif($desc_dip_taxon eq 'Bpin' || $desc_dip_taxon eq 'Bsyl'){ $lineage_code = 'G' }
#		}
#		elsif($anc_dip_taxon eq 'Bpin'){ $lineage_code = 'I' }
#		elsif($anc_dip_taxon eq 'Bsyl'){ $lineage_code = 'H' }
					
		if($lineage_code ne '-') # take all sequences with same label to choose longest afterwards
		{
			# count each code label once per taxon
			if(!defined($taxon_stats{$lineage_code})){
					$taxon_stats{$lineage_code}++;
					$taxon_stats{'all'}++;
			}

			push(@{$taxon_seqs{$lineage_code}},$node1->id());
		}

		#printf("%s %s %d %s %d %s\n",$node1->id(),$anc_dip_taxon,$anc_is_sister,$desc_dip_taxon,$desc_is_sister,$lineage_code);
	} 

	# print nr polyploid sequences to labelled MSA, max 1 per label (longest)
	foreach $lineage_code (@CODES){
		next if($lineage_code eq 'all');

		if(defined($taxon_seqs{$lineage_code})){
			foreach $header (sort {$length{$b}<=>$length{$a}} (@{$taxon_seqs{$lineage_code}})){
				
				printf(LABELMSA ">%s_%s %s\n%s\n",$taxon,$lineage_code,$header,$FASTA{$header});

				printf(LABELMSAREDUCED ">%s_%s %s\n%s\n",$taxon,$lineage_code,$header,$FASTA{$header});				

				$header2label{$header} = $lineage_code;

				last; #print "$taxon $lineage_code $header $length{$header}\n";
			}
		}
		else{
			printf(LABELMSA ">%s_%s\n%s\n",$taxon,$lineage_code,'-' x $MSAwidth);
		}
	}
	
	# print label stats of this polyploid taxon
	print "$taxon";
	foreach $lineage_code (@CODES){
		printf("\t%d",$taxon_stats{$lineage_code} || 0);
	}	print "\n";
}		

close(LABELMSAREDUCED);

close(LABELMSA);


# print Newick tree with labels for QC
my $unsorted_tree = Bio::Phylo::IO->parse(
  '-string' => $intree->as_text('newick'),
  '-format' => 'newick'
)->first();
$sorted_newick = $unsorted_tree->to_newick();
$sorted_newick =~ s/'//g;

foreach $taxon (keys(%header2label))
{
	$sorted_newick =~ s/$taxon/$taxon\_$header2label{$taxon}/;
}

open(LABELTREE,">",$labelled_tree_file);
print LABELTREE $sorted_newick;
close(LABELTREE);

# get recursively all ancestors of a node 
# returns a list of internal ids sorted from immediate to toplevel/root ancestor
sub get_all_ancestors
{
	my $node = $_[0];
	my @ancestors;

	while(defined($node->ancestor()))
	{
		$node = $node->ancestor();
	#	print "$node ".$node->internal_id()."\n"; # debug
		push(@ancestors,$node->internal_id());
	}

	#print join(',',@ancestors)."\n"; # debug

	return @ancestors;
}

# receives two sorted lists of ancestors 
# returns the MRCA (most_recent_common_ancestor)
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

# receives two sorted lists of ancestors (diploid, polyploid), a reference node (ancestor MRCA)
# and a hash reference with all known diploid MRCAs as keys
# returns a list of contiguous nodes shared between lists plus the next descendant diploid MRCA
# Only nodes younger than reference node are considered
sub get_shared_ancestors_next_diploid_MRCA
{
   my ($ref_ancestorsA, $ref_ancestorsB, $reference_node, $ref_dip_MRCA) = @_;
   my ($nodeA,$nodeB,$first_shared,%hashB);
   my @shared = ();

	# hash elements in list B (polyploid)
	# Bret 26,27,28,29
	# Bboi 21,22,23,24,25
	# Bret 20,21,22,23,24,25
	foreach $nodeB (@$ref_ancestorsB)
   {
      $hashB{$nodeB} = 1;
   }

	# iterate along A (diploid)
	# Bdis 25,26,27,28,29
	# Barb 10,11,19,20,21,22,23,24,25
	# Bpin 11,21,22,23,24,25
	$first_shared = 0; 
	foreach $nodeA (@$ref_ancestorsA)
   {
		last if($nodeA > $reference_node); # 27
		
      if($hashB{$nodeA})
      {
			# first shared descendant node of diploid, so we can later backtrack and find next MRCA
			if(scalar(@shared) == 0){ $first_shared = $nodeA }
			push(@shared, $nodeA);
      }
   }

	# loop backwards along diploid nodes starting from $first_shared 
	# until a diploid MRCA node is found
	foreach $nodeA (reverse(@$ref_ancestorsA))
   {
		next if($nodeA >= $first_shared);
		if($ref_dip_MRCA->{$nodeA}){
			unshift(@shared,$nodeA);
			last;
		}	
	}

	return @shared;
}
