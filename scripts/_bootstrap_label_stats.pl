#!/usr/bin/env perl
#
# reads all multi-tree Newick files in a folder and, for each file,
# i) extracts each individual tree
# ii) re-roots and ladderizes the individual tree
# iii) re-labels the polyploid leafs in the tree
# iv) computes overall frequencies of polyploid labels
#
# If optional argument 'rootonly' is passed steps iii and iv are skipped
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

my ($all_trees_dir,$rootonly,$taxon_label) = ('',0);

if(!$ARGV[0]){ die "# usage: $0 <folder with multi-tree Newick files> [rootonly]\n" }
else{ $all_trees_dir = $ARGV[0] }

if($ARGV[1] && $ARGV[1] eq 'rootonly'){ 
	$rootonly = 1 
}

print "# $0 $all_trees_dir $rootonly\n\n";

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

		next if($rootonly == 1);

		# iii) relabel polyploid leaf (actually we don't care about relabelled tree) 
		# and collect polyploid_taxon_label stats
		open(LABELSTATS,"$RELABELEXE -t $rooted_treefile -l|") ||
			die "# cannot run $RELABELEXE -t $rooted_treefile -l\n";

		while(my $line = <LABELSTATS>){ 
			next if($line =~ /^#/);
			if($line =~ /^(\S+)/){  
				my @stats = split(/\t/,$line);
				$taxon_label = $1;
				foreach $col (1 .. scalar(keys(%COLLABEL))){
					$boot_stats{$taxon_label}{$COLLABEL{$col}} += $stats[$col];#				}
				}
			}	
		}
		close(LABELSTATS);
		
		$count++;
	}
	close(MULTITREE);

	# iv) print this multi-tree label stats
	my ($stats_line, $n_of_data);	

	# shorten names for stats table 
	# 6L_OrthoMCL_group6293_0.7_NoGaps.fa.block.label...
	# 99985_c30921_g1_i2_chr9-30933669-30938122--.label...
	my $short_name = $multifile;
	if($multifile =~ m/^(\.*)?\.block/){
		$short_name = $1;
	}

	foreach $taxon_label (@polyconfig::polyploids_labelled){

		next if(!defined($boot_stats{$taxon_label}));

		$n_of_data = 0;
		$stats_line = '';

		$stats_line .= "$short_name\t$taxon_label";
		foreach $col (1 .. scalar(keys(%COLLABEL))){
			$stats_line .= sprintf("\t%d",		
				$boot_stats{$taxon_label}{$COLLABEL{$col}});
			if($boot_stats{$taxon_label}{$COLLABEL{$col}} > 0){
				$n_of_data++;
			} 
      } 

		# skip taxon_labels with no data
		next if($n_of_data == 0);

		# actually print some stats
		print "$stats_line\n";	
	}

	#last; # stop after 1 gene for debugging
}
