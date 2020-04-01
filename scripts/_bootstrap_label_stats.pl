#!/usr/bin/env perl
#
# reads all multi-tree Newick files in a folder and, for each file,
# i) extracts each individual tree
# ii) re-roots and ladderizes the individual tree
# iii) re-labels the polyploid leafs in the tree
# iv) computes overall frequencies of polyploid labels
#
# Input trees are usually produced by IQ-TREE and have the 
# extension .boottrees (see $EXT below)
#
# B Contreras-Moreira, R Sancho, EEAD-CSIC & EPS-UNIZAR 2018-20

use strict;
use warnings;
use Getopt::Std;
use FindBin '$Bin';
use lib "$Bin";
use polyconfig;
use polyutils;

my $EXT = '.boottrees';
my $REROOTEXE  = "$Bin/_reroot_tree.pl";
my $RELABELEXE = "$Bin/_check_lineages_polyploids.pl";

#########################################################

my ($all_trees_dir,$taxon_label);

if(!$ARGV[0]){ die "# usage: $0 <folder with multi-tree Newick files>\n" }
else{ $all_trees_dir = $ARGV[0] }

# read all input files names
opendir(ALL,$all_trees_dir) || die "# cannot list $all_trees_dir\n";
my @multiNewick = grep{/$EXT/} readdir(ALL);
closedir(ALL);

# hash positions of @polyconfig::CODES
my (%COLLABEL, $col);
foreach $col (0 .. $#polyconfig::CODES) {
   $COLLABEL{ $col+1 } = $polyconfig::CODES[ $col ];
}

# parse multi-tree files one at a time
foreach my $multifile (sort (@multiNewick)){
	
	my (%boot_stats);	
	
	# create folder to save individual trees 
	my $resultsdir = "$all_trees_dir/$multifile.relabelled";
	mkdir($resultsdir);
	if(!-e $resultsdir){ die "# cannot create $resultsdir\n"; }

	# parse individual trees
	my $count = 1;
	open(MULTITREE,"$all_trees_dir/$multifile") || 
		die "# cannot read $all_trees_dir/$multifile\n";

	while(<MULTITREE>){
		next if(/^$/);
		
		# i) extract individual tree to file
		my $treefile = "$resultsdir/$count.ph";	
		open(PH,">",$treefile) || die "# cannot create $treefile\n";
		print PH;
		close(PH);

		# ii) re-root and ladderize individual trees
		my $rooted_treefile = "$resultsdir/$count.root.ph";
		system("$REROOTEXE $treefile > $rooted_treefile");

		# iii) relabel polyploid leaf (actually we don't care about relabelled tree) 
		# and collect polyploid_taxon_label stats
		print "$RELABELEXE -t $rooted_treefile\n"; exit;
		open(LABELSTATS,"$RELABELEXE -t $rooted_treefile |") ||
			die "# cannot run $RELABELEXE -t $rooted_treefile\n";

		while(<LABELSTATS>){ 
			if(/^>(\S+)/){ print;
				my @stats = split(/\t/,$_);
				$taxon_label = $1;
				foreach $col (1 .. scalar(keys(%COLLABEL))){
					$boot_stats{$taxon_label}{$COLLABEL{$col}} += $stats[$col];
				}
			}	
		}
		close(LABELSTATS);
		exit;
		$count++;
	}
	close(MULTITREE);

	# iv) print this multi-tree label stats
	# 6L_OrthoMCL_group6293_0.7_NoGaps.fa.block.label.reduced.fna.trimmed.fna.extract_Ttur_I.fna.boottrees
	my $short_name = $multifile;
	if($multifile =~ m/^(\.*_group\d+)/){
		$short_name = $1;
	}
	foreach $taxon_label (@polyconfig::polyploids_labelled){

		next if(!defined($boot_stats{$taxon_label}));

		print "$short_name\t$taxon_label";
		foreach my $col (1 .. 9){
			printf("\t%d",		
				$boot_stats{$taxon_label}{$COLLABEL{$col}});
      } print "\n";
	}

	#last; # stop after 1 gene
}
