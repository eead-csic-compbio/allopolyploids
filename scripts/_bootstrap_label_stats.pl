#!/usr/bin/env perl
use strict;
use warnings;

# reads all multi-tree Newick files in a folder and, for each file,
# i) extracts each individual tree
# ii) re-roots and ladderizes the individual tree
# iii) re-labels the polyploid leafs in the tree
# iv) computes overall frequencies of polyploid labels

my %COLLABEL = ( 1,'A',2,'B',3,'C',4,'D',5,'E',6,'F',7,'G',8,'H',9,'I' );
my @polyploids = (
        'Bhyb_A','Bhyb_B','Bhyb_C','Bhyb_D','Bhyb_E','Bhyb_F','Bhyb_G','Bhyb_H','Bhyb_I',
        'Bboi_A','Bboi_B','Bboi_C','Bboi_D','Bboi_E','Bboi_F','Bboi_G','Bboi_H','Bboi_I',
        'Bret_A','Bret_B','Bret_C','Bret_D','Bret_E','Bret_F','Bret_G','Bret_H','Bret_I',
        'Bmex_A','Bmex_B','Bmex_C','Bmex_D','Bmex_E','Bmex_F','Bmex_G','Bmex_H','Bmex_I',
        'Brup_A','Brup_B','Brup_C','Brup_D','Brup_E','Brup_F','Brup_G','Brup_H','Brup_I',
        'Bpho_A','Bpho_B','Bpho_C','Bpho_D','Bpho_E','Bpho_F','Bpho_G','Bpho_H','Bpho_I',
        'B422_A','B422_B','B422_C','B422_D','B422_E','B422_F','B422_G','B422_H','B422_I'
);


my $EXT = '.boottrees';
my $REROOTEXE = './_reroot_tree.pl';
my $RELABELEXE = './_check_lineages_polyploids_ABCDEFGHI_noFASTA.pl';

my ($all_trees_dir,$taxon_label);

if(!$ARGV[0]){ die "# usage: $0 <folder with multi-tree Newick files>\n" }
else{ $all_trees_dir = $ARGV[0] }

opendir(ALL,$all_trees_dir) || die "# cannot list $all_trees_dir\n"; 
my @multiNewick = grep{/$EXT/} readdir(ALL);
closedir(ALL);

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
		open(LABELSTATS,"$RELABELEXE $rooted_treefile |") ||
			die "# cannot run $RELABELEXE $rooted_treefile\n";
		while(<LABELSTATS>){
			if(/^>(\S+)/){ #print;
				my @stats = split(/\t/,$_);
				$taxon_label = $1;
				foreach my $col (1 .. 9){
					$boot_stats{$taxon_label}{$COLLABEL{$col}} += $stats[$col];
				}
			}	
		}
		close(LABELSTATS);

		$count++;
	}
	close(MULTITREE);

	# iv) print this multi-tree label stats
	my $short_name = $multifile;
	if($multifile =~ m/^(\d+_c\d+)/){
		$short_name = $1;
	}
	foreach $taxon_label (@polyploids){

		next if(!defined($boot_stats{$taxon_label}));

		print "$short_name\t$taxon_label";
		foreach my $col (1 .. 9){
			printf("\t%d",		
				$boot_stats{$taxon_label}{$COLLABEL{$col}});
      } print "\n";
	}

	#last; # stop after 1 gene
}
