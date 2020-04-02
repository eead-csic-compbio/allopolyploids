#!/usr/bin/perl 

use strict;
use warnings;
use Getopt::Std;
use FindBin '$Bin';
use lib "$Bin";
use polyconfig;
use polyutils;

# external dependencies, need to be installed
use IO::String;
use Bio::TreeIO;
use Bio::Phylo::IO;
use Bio::Phylo::Forest::Tree;

# Label allopolyploid sequences with pre-defined codes according to their position
# with respect to the diploid backbone. When several sequences of the same
# allopolyploid species are found, all of them are labelled.
# Please see and edit polyconfig.pm to define species, codes and rules. 
#  
# INPUT: i) multiple alignment (MSA, .fna) corresponding to ii) tree (.treefile.root.ph) 
# OUTPUT: produces three output files with labelled sequences:
# i) full MSA with sorted diploid/polyploid species, with gap-only rows, can be concatenated
# ii) reduced MSA excluding gap-only rows
# iii) tree in Newick format
#
# Ad-hoc rules encoded in module polyconfig.pm define how labels are assigned based on
# the ascendant and descendant diploid species.
#
# Optionally sister clades might be defined and their most recent common ancestors (MRCA)
# used in the labelling rules (see polyconfig.pm).
#
# B Contreras-Moreira, R Sancho, EEAD-CSIC & EPS-UNIZAR 2018-20

my ($outfound,$outnode,$sorted_newick,$taxon,$node,$node_id,$child,$anc) = (0);
my ($min_anc_dist,$min_desc_dist,$dist,$anc_dip_taxon,$desc_dip_taxon) = (0,0);
my $MSAwidth = 0;
my ($coord,$coord5,$coord3,$midcoord);
my ($d,$header,$node1,$node2,$taxon1,$taxon2,$bases,$lineage_code,$full_taxon);
my ($poly_ancestor_MRCA,$desc_is_sister,$anc_is_sister,$next_MRCA);
my (%length,%FASTA,%order,%dip_taxon,%MRCA_computed,%dip_ancestor,%header2label);
my (%dip_all_ancestors,%dip_all_ancestors_string,%diptaxon2node);
my ($MRCA,$diploid_MRCA,$diploid_MRCA_node,%dip_MRCA);
my (@clade_MRCA_nodes,%clade_MRCA,%clade_ancestors);

my ($MSA_file,$tree_file,$use_labeled_polyploids,$verbose,%opts) = ('','',0,0);

getopts('hf:t:lv', \%opts);

if(($opts{'h'})||(scalar(keys(%opts))==0))
{
	print "\nusage: $0 [options]\n\n";
	print "-h this message\n";
	print "-t rooted, sorted tree Newick file with .ph extension\n";   
	print "-f FASTA file with aligned tree sequences and .fna extension (optional)\n";
   print "-v verbose output, useful for debugging                      (optional)\n";
	print "-l use CODE-labelled polyploids instead of plain names,      (optional)\n";
	print "   this is useful to relabel boostrapped trees\n\n";
	exit(0);
} 

if(defined($opts{'t'}) && -s $opts{'t'} && $opts{'t'} =~ /\.ph/){ 
	$tree_file = $opts{'t'}; 
}
else{ die "# EXIT : need a valid tree file with .ph extension\n"; }

if(defined($opts{'f'}) && -s $opts{'f'} && $opts{'f'} =~ /\.fna/){
   $MSA_file = $opts{'f'};
}

if(defined($opts{'v'})){
   $verbose = 1;
}

if(defined($opts{'l'})){
	$use_labeled_polyploids = 1;
}

# print parsed arguments
print "# $0 -t $tree_file -f $MSA_file -v $verbose -l $use_labeled_polyploids\n\n";

# set polyploid list of names to use
my @polyploids = @polyconfig::polyploids;
if($use_labeled_polyploids){
	@polyploids = @polyconfig::polyploids_labelled;
}

# names of three output files
my $labelledMSA_file = $MSA_file;
$labelledMSA_file =~ s/\.fna/.label.fna/;
my $reducedMSA_file = $MSA_file;
$reducedMSA_file =~ s/\.fna/.label.reduced.fna/;
my $labelled_tree_file = $tree_file;
$labelled_tree_file =~ s/\.ph/.label.ph/;

# read FASTA sequences and record length
if($MSA_file ne ''){
	open(MSA,"<",$MSA_file) || die "# ERROR: cannot read $MSA_file\n"; 
	while(<MSA>){
		if(/^>(\S+)/){ 
			$header = $1; 
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
}

# ladderize/sort input tree 
my $unsorted_input_tree = Bio::Phylo::IO->parse(
  '-file' => $tree_file,
  '-format' => 'newick'
)->first();

$unsorted_input_tree->ladderize($polyconfig::NODEORDER);
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
my (@refDiploid, %dip_nodes);
my $total_nodes = 0;
for $node ($intree->get_nodes()) {
	if(defined($node->id())) {

		$total_nodes++; # labelled nodes that is

		foreach $d ( 0 .. $#polyconfig::diploids) {			
			if($node->id() =~ m/_($polyconfig::diploids[$d])$/ || # normal trees
					$node->id() =~ m/^($polyconfig::diploids[$d])$/) { #label bootstrapped

				$taxon1 = $1;
				push(@refDiploid,$node); 
				$dip_taxon{$node} = $taxon1;
				$diptaxon2node{$taxon1} = $node; 
				$order{$node->id()} = $total_nodes; 
				$dip_ancestor{$node} = $node->ancestor(); # immediate ancestor
				print "$total_nodes ".$node->id()." ".$node->ancestor()."\n" if($verbose);

				my @all_ancestors = get_all_ancestors($node);
				foreach $anc (@all_ancestors){ $dip_nodes{$anc}++ }				
				$dip_all_ancestors_string{$node} = join(',',@all_ancestors);
				$dip_all_ancestors{$node} = \@all_ancestors;  
				print "$dip_all_ancestors_string{$node}\n" if($verbose);
			}
		}
	}	
}

if(scalar(@refDiploid) == 0){
	die "# ERROR: failed to parsed any diploids, please check their names\n"
} 

# identify most recent common ancestors (MRCA) among diploid nodes
# diploids nodes are compared pairwise in the order they appear from root to tips
foreach $node (0 .. $#refDiploid-1){

	$node1 = $refDiploid[$node];
	$node2 = $refDiploid[$node+1]; 

   if(!defined($MRCA_computed{$node1})){ # sister of previous node

		# assumes node 1 is ancestor to node 2, 
		# note that $MRCA is an internal node id, which is an arbitrary integer
		$MRCA = get_MRCA( $dip_all_ancestors{$node1} , $dip_all_ancestors{$node2} ); 
		$dip_MRCA{ $MRCA } = $dip_taxon{$node1};
		$MRCA_computed{$node1} = 1;
		print "MRCA $dip_MRCA{ $MRCA } $MRCA\n" if($verbose); 

		# but we should check whether they are sisters
		if($dip_all_ancestors_string{$node1} eq $dip_all_ancestors_string{$node2}){
			print "MRCA $dip_taxon{ $node2 } $MRCA (sister)\n" if($verbose);
			$MRCA_computed{$node2} = 1;
		}	
   }

	# add also ancestor of last diploid in order
	if($node == $#refDiploid-1 && !defined($MRCA_computed{$node2})){
		$MRCA = $node2->ancestor()->internal_id(); # fake MRCA for last diploid
		$dip_MRCA{ $MRCA } = $dip_taxon{$node2};
      $MRCA_computed{$node2} = 1;
		print "MRCA $dip_MRCA{ $MRCA } $MRCA\n" if($verbose); 
	}
}


## search for MRCA nodes of any sister clades defined above
## and save them in %clade_MRCA and %clade_ancestors
my ($snode1,$snode2,$snode,$sclade,$sMRCA,$MRCAfound,%sisterMRCA);
foreach my $bifur (keys(%polyconfig::sister_clades)){

	my %overwrite;

	# find MRCA of each clade and save it
	foreach $sclade (keys(%{ $polyconfig::sister_clades{$bifur} })){
		
		my (@sancestors);
		my @snodes = @{ $polyconfig::sister_clades{$bifur}->{$sclade} };
		$sMRCA = -1; # init

		foreach $snode (0 .. $#snodes-1){
			$taxon1 = $snodes[$snode];
			$taxon2 = $snodes[$snode+1]; #print "$taxon1 $taxon2\n";
			
			# 1st MRCA is computed for this clade
			if($sMRCA == -1){ 
				$sMRCA = get_MRCA( $dip_all_ancestors{$diptaxon2node{$taxon1}} , 
								$dip_all_ancestors{$diptaxon2node{$taxon2}} );
			
				# retrieve the list of node ids of the ancestors, 
				# which should be among the ancestors of taxon2
				$MRCAfound = 0;
				foreach $anc ( @{ $dip_all_ancestors{ $diptaxon2node{$taxon2} } } ){
					if($anc == $sMRCA){ $MRCAfound = 1 } # ensure MRCA is added as ancestor
					if($MRCAfound == 1){ push(@sancestors,$anc) }
				} 
			}
			else { #subsequent MRCA calculations within clade
				$sMRCA = get_MRCA( \@sancestors , 
								$dip_all_ancestors{$diptaxon2node{$taxon2}} );

				# retrieve again the list of node ids of the ancestors,
				@sancestors = ();
				$MRCAfound = 0;
            foreach $anc ( @{ $dip_all_ancestors{ $diptaxon2node{$taxon2} } } ){
               if($anc == $sMRCA){ $MRCAfound = 1 }
					if($MRCAfound == 1){ push(@sancestors,$anc) }
            }
			}
		}
		
		# save this clade's MRCA (node) & ancestors
		if(defined($clade_MRCA{$sMRCA})){
			print "# overwrite $clade_MRCA{$sMRCA} with $sclade\n" if($verbose);
			$overwrite{$sclade}{$clade_MRCA{$sMRCA}} = 1; # Tura,Tmon
			$overwrite{$clade_MRCA{$sMRCA}}{$sclade} = 1; # Tmon,Tura
		}
		$clade_MRCA{$sMRCA} = $sclade;
		$clade_ancestors{$sMRCA} = \@sancestors;
		print "MRCA $sclade $sMRCA\n" if($verbose);
	}

	# find MRCA of the two clades in this bifurcation,
	# note that MRCA may well be one of those already computed
	my @sclades = keys(%{ $polyconfig::sister_clades{$bifur} });
	($node1,$node2) = ('','');
	foreach $snode (keys(%clade_MRCA)){
		if($clade_MRCA{$snode} eq $sclades[0]){ $node1 = $snode }
		elsif($clade_MRCA{$snode} eq $sclades[1]){ $node2 = $snode }
	}

	# in case MRCA nodes 1/2 were overwritten
	if(!$node1){ 
		if($overwrite{$sclades[0]}{$sclades[1]}){ $node1 = $node2 } 
		else{ die "# ERROR: cannot find node $sclades[0]\n" }
	}
   elsif(!$node2){ 
		if($overwrite{$sclades[0]}{$sclades[1]}){ $node2 = $node1 } 
		else{ die "# ERROR: cannot find node $sclades[1]\n" }
	}

	$sMRCA = get_MRCA( $clade_ancestors{ $node1 } ,
                         $clade_ancestors{ $node2 } );

	#print "<$bifur> $sMRCA $node1 $node2 ".
   #               join(',',@{$clade_ancestors{$node1}})." ".
   #               join(',',@{$clade_ancestors{$node2}})."\n" if($verbose);

	# find its ancestors
	$MRCAfound = 0;
	my @cancestors;
   foreach $anc ( @{ $clade_ancestors{ $node1 } } ){
		if($MRCAfound == 1){
			push(@cancestors,$anc);
      }
      if($anc == $sMRCA){ $MRCAfound = 1 }
   }

	# save this sister clade MRCA (node) & ancestors
	if(defined($clade_MRCA{$sMRCA})){
      print "# overwrite $clade_MRCA{$sMRCA} with $bifur\n" if($verbose);
   }
	$clade_MRCA{$sMRCA} = $bifur;
	$clade_ancestors{$sMRCA} = \@cancestors;
	print "MRCA $bifur $sMRCA\n" if($verbose);

	# save list of nodes in hierarchical order
	push(@clade_MRCA_nodes,$sMRCA);
	if($node1 != $sMRCA){ push(@clade_MRCA_nodes,$node1) } 
	if($node2 != $sMRCA){ push(@clade_MRCA_nodes,$node2) }
} 

# create labelled MSA files
if($MSA_file ne ''){
	open(LABELMSA,">",$labelledMSA_file) || 
		die "# ERROR: cannot create $labelledMSA_file\n";

	open(LABELMSAREDUCED,">",$reducedMSA_file) || 
		die "# ERROR: cannot create $reducedMSA_file\n";

	# write output file names
	print "\n# $labelledMSA_file $reducedMSA_file $MSAwidth $labelled_tree_file\n\n";
} else {
	print "\n# NA NA NA $labelled_tree_file\n\n";
}


# print diploid sequences
if($MSA_file ne ''){
	foreach $taxon (@polyconfig::diploids){
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
		if(!$matched) {
			printf(LABELMSA ">%s\n%s\n",$taxon,'-' x $MSAwidth);
		}
	}
}

# print header of table of labels (see @CODES)
foreach $lineage_code (@CODES){ printf("\t%s",$lineage_code) }  print "\n";

# find polyploid nodes, label nr and add them to labelled MSA
foreach $taxon (@polyploids)
{
	my (%taxon_ancestors,@nr_taxon_nodes,%taxon_order,%taxon_stats,%taxon_seqs);

	print "> $taxon\n" if($verbose);

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

	# simplify same-taxon polytomies
	foreach $anc (keys(%taxon_ancestors)){
		my $total_polytomic_nodes = 0;
		my @polytomic_nodes;

		# if sequences are known, take longest sequence
		if($MSA_file ne ''){
			@polytomic_nodes = 
				sort {$length{$b->id()}<=>$length{$a->id()}} @{$taxon_ancestors{$anc}};			
		}
		else { # if sequences are not known, choose 1st one
			@polytomic_nodes = @{$taxon_ancestors{$anc}};
      }

		foreach $node (@polytomic_nodes){
			$total_polytomic_nodes++;
			if($total_polytomic_nodes == 1){
				push(@nr_taxon_nodes,$node); 
				#printf("%s %d\n",$node->id(),$length{$node->id()});
			}
		} 
	}

	# check closest diploid ancestor and descendant of this polyploid node
	foreach $node1 (@nr_taxon_nodes){

		($anc_dip_taxon,$desc_dip_taxon,$lineage_code) = ('-','-','-');
		$poly_ancestor_MRCA = -999;
		$anc_is_sister = 0;

		# get all ancestors of this polyploid node, from immediate to toplevel/root
		my @node1_all_ancestors = get_all_ancestors($node1); 
		my $node1_all_ancestors_string = join(',',@node1_all_ancestors);
		print $node1->id()." $node1_all_ancestors_string\n" if($verbose); 

		# find closest ancestor which is MRCA of two diploids (or fake MRCA of last diploid)
		my $total_dip_nodes = 0;
		foreach $node_id (@node1_all_ancestors){
	
			# check whether this node is 1st ancestor of a diploid node,
			# meaning it is sister to this diploid
			if($dip_nodes{$node_id}){
				$total_dip_nodes++;
				foreach $node2 (@refDiploid){
					if($total_dip_nodes ==1 && $node_id == $dip_all_ancestors{$node2}->[0]){
						$poly_ancestor_MRCA = $node_id;
						$anc_dip_taxon = $dip_taxon{ $node2 };
						$anc_is_sister = 1;
						last;
					}
				}
				last if($anc_is_sister == 1);
			}

			#otherwise check whether this node is a previously defined clade/diploid MRCA
			if(defined($clade_MRCA{ $node_id })){
            $poly_ancestor_MRCA = $node_id;
            $anc_dip_taxon = $clade_MRCA{ $node_id };
				last;
         }
			elsif(defined($dip_MRCA{ $node_id })){
				$poly_ancestor_MRCA = $node_id;
				$anc_dip_taxon = $dip_MRCA{ $node_id }; 
				last;
			}
		}
		print "anc_dip_taxon $anc_dip_taxon ($anc_is_sister)\n" if($verbose);		
	
		# find the first descendant which is either: 
		# i) a MRCA of sister clades, ii) a diploid sister, iii) a MRCA between diploids
		my (%total_shared_desc,%desc_sister, @sorted_taxa); 

		# try first MRCA nodes of explicitely defined sister nodes
		foreach $node2 (@clade_MRCA_nodes){

			# a taxon cannot be both ancestor and descendant
			next if($anc_dip_taxon eq $clade_MRCA{$node2});

			#$next_MRCA = '';
			my @shared_nodes = get_shared_ancestors( $clade_ancestors{$node2},
					\@node1_all_ancestors, $poly_ancestor_MRCA );

			#if(scalar(@shared_nodes) > 1){		 
			#	$next_MRCA = get_next_MRCA( $clade_ancestors{$node2}, $shared_nodes[0], 
			#		\%clade_MRCA  ); 
			#} 

			print "<$clade_MRCA{$node2}> ".join(',',@shared_nodes).
                  " ($poly_ancestor_MRCA) ".
                  join(',',@{$clade_ancestors{$node2}})." ".
                  join(',',@node1_all_ancestors)."\n" if($verbose);	

			$desc_dip_taxon = $clade_MRCA{$node2};
         $total_shared_desc{$desc_dip_taxon} = scalar(@shared_nodes);
			push(@sorted_taxa, $clade_MRCA{$node2});
		} 
	
		# now try diploids
		foreach $node2 (@refDiploid){
		
			# a taxon cannot be both ancestor and descendant
         next if($anc_dip_taxon eq $dip_taxon{$node2});

			# all ancestors are shared: polyploid (node1) and diploid (node2) are sisters
			if($node1_all_ancestors_string eq $dip_all_ancestors_string{$node2}){

				$desc_dip_taxon = $dip_taxon{ $node2 };
				$total_shared_desc{$desc_dip_taxon} = scalar(@{$dip_all_ancestors{$node2}});
				$desc_sister{$desc_dip_taxon} = 1;
				push(@sorted_taxa, $desc_dip_taxon);
				last; #print "des_tip_taxon $desc_dip_taxon\n" if($verbose); 
			}
			else{ # check shared ancestors between consecutive diploid MRCA nodes

				$next_MRCA='';
				my @shared_nodes = get_shared_ancestors( $dip_all_ancestors{$node2}, 
								\@node1_all_ancestors, $poly_ancestor_MRCA );

				if(scalar(@shared_nodes) > 1){
					$next_MRCA = get_next_MRCA( $dip_all_ancestors{$node2}, $shared_nodes[0],
						\%dip_MRCA  );
				}

				print "<$dip_taxon{ $node2 }> ".join(',',@shared_nodes).
						" ($poly_ancestor_MRCA) ". 
						join(',',@{$dip_all_ancestors{$node2}})." ".
						join(',',@node1_all_ancestors)."\n" if($verbose);

				if(scalar(@shared_nodes) > 0){
								
					if($next_MRCA ne ''){ 
						$desc_dip_taxon = $dip_MRCA{$next_MRCA};
						$total_shared_desc{$desc_dip_taxon} = scalar(@shared_nodes);
						#print "des_dip_taxon $desc_dip_taxon $shared_nodes[0]\n" if($verbose);
					}
					else{ # in case there are no more diploid MRCAs down the line
		            $desc_dip_taxon = $dip_taxon{ $node2 };
						$total_shared_desc{$desc_dip_taxon} = scalar(@shared_nodes);
						$desc_sister{$desc_dip_taxon} = 1;
				      #print "des_tip_taxon $desc_dip_taxon\n" if($verbose);
					}

					push(@sorted_taxa, $desc_dip_taxon);
				}
			}
		}

		# check diploid with most shared ancestors, which will be descendant
		my $max_shared = 0;
		$desc_is_sister = 0;	
		$desc_dip_taxon = '';
		foreach $taxon (@sorted_taxa){
			if(defined($total_shared_desc{$taxon}) && 
					$total_shared_desc{$taxon} > $max_shared){
				$desc_dip_taxon = $taxon;
				$max_shared = $total_shared_desc{$taxon};
			} 
		}

		# is descendant a sister node?
		if(defined($desc_sister{$desc_dip_taxon}) && 
			$desc_sister{$desc_dip_taxon} == 1){
			$desc_is_sister = 1;
		}
	
		$lineage_code = get_label_from_rules( $anc_dip_taxon, $desc_dip_taxon, 
															$anc_is_sister, $desc_is_sister);

		print "des_dip_taxon $desc_dip_taxon ($desc_is_sister) -> $lineage_code\n" if($verbose);

		# collect some label stats
		if($lineage_code ne '-')
		{
			# count each code label once per taxon
			if(!defined($taxon_stats{$lineage_code})){
					$taxon_stats{$lineage_code}++;
					$taxon_stats{'all'}++;
			}

			push(@{$taxon_seqs{$lineage_code}},$node1->id());
		}

	} 

	# add labels to sequence headers and tree nods
	# print nr polyploid sequences to labelled MSA, max 1 per label
	foreach $lineage_code (@CODES){
		next if($lineage_code eq 'all');

		if(defined($taxon_seqs{$lineage_code})){

			my @headers = @{$taxon_seqs{$lineage_code}};
			if($MSA_file ne ''){	
				@headers = sort {$length{$b}<=>$length{$a}} (@{$taxon_seqs{$lineage_code}});
			}

			foreach $header (@headers){

				if($MSA_file ne ''){	
					printf(LABELMSA ">%s_%s %s\n%s\n",
						$taxon,$lineage_code,$header,$FASTA{$header});

					printf(LABELMSAREDUCED ">%s_%s %s\n%s\n",
						$taxon,$lineage_code,$header,$FASTA{$header});				
				}
				
				$header2label{$header} = $lineage_code;

				# one first sequence considered per taxon
				# if sequence lengths are available, it would be the longest one
				last; 
			}
		}
		else{
			if($MSA_file ne ''){
				printf(LABELMSA ">%s_%s\n%s\n",$taxon,$lineage_code,'-' x $MSAwidth);
			}
		}
	}

	
	# print label stats of this polyploid taxon
	print "$taxon";
	foreach $lineage_code (@CODES){
		printf("\t%d",$taxon_stats{$lineage_code} || 0);
	}	print "\n";
}		

if($MSA_file ne ''){
	close(LABELMSAREDUCED);

	close(LABELMSA);
}


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
