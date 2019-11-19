# allopolyploids

This pipeline was designed by Ruben Sancho, Antonio Diaz, Pilar Catalan and Bruno Contreras Moreira for the selection of transcripts for phylogeny reconstruction of allopolyploid species. It uses what we call the **nearest diploid species node algorithm**. We tested it with diploid and polyploid species of the genus Brachypodium, for which we had data obtained in collaboration with David des Marais. The ideas and the code can be taylored for other clades as well.


## 1) Core transcripts expressed in all Brachypodium species plus two outgroups: rice & barley

This first step requires https://github.com/eead-csic-compbio/get_homologues and the set of transcripts in folder [genome_transcripts](./genome_transcripts), most of them assembled de novo with https://github.com/trinityrnaseq/trinityrnaseq

```
get_homologues/get_homologues-est.pl -d genome_transcripts/ -m cluster -I \
	genome_transcripts/species.list -M -A -S 80 &> \
	00_get_homologues/log.gen.M.A.S80.core.clusters
```

This produces 3324 clusters (see this [folder](./00_get_homologues/genome_transcripts_est_homologues/arb8075_alltaxa_species.list_algOMCL_e0_S80_)) and an Average Nucleotide Identity (ANI) matrix, which we can plot with: 

```
get_homologues/plot_matrix_heatmap.sh -i 00_get_homologues/genome_transcripts_est_homologues/arb8075_alltaxa_species.list_algOMCL_e0_S80_Avg_identity.tab \
	-H 10 -W 18 -t "ANI of transcripts in 3324 core clusters" -o svg -d 1
```
![ANI matrix](00_get_homologues/genome_transcripts_est_homologues/arb8075_alltaxa_species.list_algOMCL_e0_S80_Avg_identity_heatmap.svg)

We can now compute pan-gene matrices for these core clusters:
```
get_homologues/compare_clusters.pl -d 00_get_homologues_est/genome_transcripts_est_homologues/arb8075_alltaxa_species.list_algOMCL_e0_S80_ \
   -o 00_get_homologues/core_clusters_Hordeum -m -n &> 00_get_homologues/log.compare.core
```

**Note**: this creates folder **core_cluster_Hordeum**

## 2) Multiple alignment trimmed blocks of core clusters

From the clusters obtained earlier we can now produce [MVIEW](https://github.com/desmid/mview) collapsed, multiple sequence alignments (MSA), while also annotating Pfam domains:
```
for FILE in `ls 00_get_homologues/core_clusters_Hordeum/*gethoms.fna; do
   echo $FILE;
   get_homologues/annotate_cluster.pl -D -f $FILE -o $FILE.aln.fna -c 20 &>> \
		00_get_homologues/log.core.collapse.align
done
```

We might need to simplify the FASTA headers, which is something that should be taylored according to the user's data. The resulting files are stored in folder **01_core_transcript_aligned**: 
```
for FILE in *collapsed.fna.gethoms.fna.aln.fna; do
echo $FILE;
perl -p -i -e 's/>(.+?) .+/>$1/g; s/:\d+:\d+:[+-]//g' $FILE;
done
```

We now take these alignments of nucleotide sequences of both diploid and polyploid species and produce trimmed FASTA files suitable for phylogenetic tree inference. The goal is to define a solid diploid backbone, which should be covered by outgroup sequences as well, and then use it to filter out polyploid sequences/alleles with diploid block overlap < $MINBLOCKOVERLAP. Therefore, some sequences and species might be removed from the initial input. These parameters are set in [scripts/_trim_MSA_block.pl](./scripts/_trim_MSA_block.pl), which also shortens the species names:
```
my $MINBLOCKLENGTH = 100;
my $MINBLOCKOVERLAP = 0.50; # fraction of diploid block covered by outgroups and polyploid seqs

# short names of species used to define diploid block (see %long2short below)
# 1 per taxon will be selected for optimizing diploid block
my @diploids = qw( _Bsta _Bdis _Barb _Bpin _Bsyl );

# diploid outgroups: best block-overlapping sequence will be conserved
my @outgroups = qw( _Osat _Hvul );
```
The resulting files are stored in folder **02_blocks**:
```
for FILE in *.fna; do
echo $FILE;
perl scripts/_trim_MSA_block.pl $FILE $FILE.block.fna &>> log.blocks;
done
```

We can now trim the resulting blocks with https://vicfero.github.io/trimal . The results are stored in folder [03_blocks_trimmed](./03_blocks_trimmed). **Note** that the number valid MSA has now reduced to just over 1700:
```
for FILE in *block.fna; do
echo $FILE;
trimal/source/trimal -in $FILE -out $FILE.trimmed.fna -automated1;
done
```

## 3) Maximum Likelihood (ML) gene trees 

We can now run [IQ-TREE](http://www.iqtree.org) in [parallel](https://www.gnu.org/software/parallel) 
and store the results in folder **04_iqtree**:
```
ls *.trimmed.fna | parallel --gnu -j 3 iqtree-omp-1.5.5-Linux/bin/iqtree-omp -alrt 1000 -bb 1000 -nt 3 -AICc -s {} :::
```

We will now root and ladderize the nodes (starting with rice) in the resulting trees, which are stored in folder [05_iqtree_rooted_sorted](./05_iqtree_rooted_sorted). **Note** this script requires Perl modules [Bio::TreeIO](https://metacpan.org/pod/Bio::TreeIO) and 
[Bio::Phylo](https://metacpan.org/pod/Bio::Phylo):
```
sudo cpan -i Bio::TreeIO Bio::Phylo

for FILE in *treefile; do
perl _reroot_tree.pl $FILE > $FILE.root.ph;
echo $FILE;
done 
```

We can now ask how many different diploid backbones are there and select gene trees with consistent diploid backbone. The next script must be run within folder **05_iqtree_rooted_sorted** and requires http://cegg.unige.ch/newick_utils :
```
for FILE in *root.ph; do
perl _check_diploids.pl $FILE _Osat; done > log.diploids;
done

grep -v "#" log.diploids | \
   perl -F"," -ane 'foreach $tx (0 .. $#F){ $ord{$F[$tx]}{$tx}++ } END{ foreach $tx (keys(%ord)){ print "$tx"; foreach $t (0 .. 6){ printf("\t%d",$ord{$tx}{$t}||0) } print "\n" } }'
```

We obtain these statistics (**note:** * highlights topologies with sister diploids):
```
Osat    1707    0       0       0       0       0       0
Hvul    0       1527    46      31      17      7       0
Bsta    0       64      936     247     90      53      0
Bdis    0       54      444     708     145     124     0
Barb    0       22      75      104     585     596     0
Bpin    0       20      72      100     390     567     167
Bsyl    0       20      77      112     363     286     431

Hvul*   0       0       10      11      0       0       18
Bsta*   0       0       15      198     20      9       75
Bdis*   0       0       12      128     23      7       62
Barb*   0       0       5       19      30      26      245
Bpin*   0       0       8       25      24      8       294
Bsyl*   0       0       7       24      20      23      332
```

And consequently save our topology predominance list in file [list.diploids_congruent_pruned_diploid_topology](./05_iqtree_rooted_sorted/list.diploids_congruent_pruned_diploid_topology):
```
sort log.diploids  | uniq -c | sort -n

# NOTE: not all grep versions have --no-group-separator flag
grep -B 8 --no-group-separator -e 'Osat,Hvul,Bsta,Bdis,Barb,Bsyl,Bpin,' -e 'Osat,Hvul,Bsta,Bdis,Barb,Bpin,Bsyl,' log.diploids | grep treefile.root.ph | sed 's/# //g' | sed 's/ /\n/g' > list.diploids_congruent_pruned_diploid_topology
```

We'll now copy files (.fna, root.ph, root.ph.pruned) with congruent diploid topologies to folder [06_diploid_clusters_pruned_diploid_topology](./06_diploid_clusters_pruned_diploid_topology):
```
ls | grep -f list.diploids_congruent_pruned_diploid_topology | \
	xargs cp -t 06_diploid_clusters_pruned_diploid_topology/

ls *.root.ph | sed 's/.treefile.root.ph//g' > list_fna_06_diploid_clusters_pruned_diploid_topology

ls | grep -f list_fna_06_diploid_clusters_pruned_diploid_topology | \
	xargs cp -t 06_diploid_clusters_pruned_diploid_topology/
```

## 4) Trees and MSA with topology-labelled allopolyploid sequences

In this step we'll label allopolyploid sequences with some pre-defined code according to their position with respect to the diploid backbone. We'll do that with [scripts/_check_lineages_polyploids_ABCDEFGHI.pl](./scripts/_check_lineages_polyploids_ABCDEFGHI.pl), which has several parameters defined therein:

```
my @diploids = qw( Osat Hvul Bsta Bdis Barb Bpin Bsyl ); # expects one seq per species and tree
my @polyploids = ('Bhyb','Bboi','Bret','Bmex','Brup','Bpho','B422');
my @CODES = qw( A B C D E F G H I all ); # see hard-coded rules below
```
This script takes as input an MSA corresponding to an input tree file and produces two labelled MSA files with sorted diploid and polyploid (A,B,C,etc) species, a full MSA with some only-gap rows, and a reduced MSA excluding those. Ad-hoc rules encoded herein (see below) define pairs of diploid species where one on the can be missing in the alignment. These rules should be carefully adapted to other clades. Labelled files are stored in folder [07_files_labelled_ABCDEFGHI](./07_files_labelled_ABCDEFGHI): 

```

ls -1 *root.ph | perl -lne 'print `../../_check_lineages_polyploids_ABCDEFGHI.pl $_`' &> log.lineage_codes_blocks_ABCDEFGHI
````

**Important:** labelled files must be curated manually and any errors corrected

In our study we obtained these label stats:

```
IMPORTANT: these labels are in lowercase in the published paper

      A    B     C    D    E    F    G    H    I     all
Bhyb  2    137   2    90   0    0    0    0    1     232
Bmex  94   46    42   4    15   1    2    0    0     204
Bboi  87   37    49   9    32   5    4    0    0     223
Bret  42   12    40   9    63   21   31   9    17    244
Brup  0    1     2    5    53   29   71   49   45    255
B422  1    2     3    4    56   40   66   37   36    245
Bpho  1    0     4    4    66   42   64   30   43    254
```

## 5) Validated labelled trees and MSA files 

Starting with the input files in [07_files_labelled_ABCDEFGHI](./07_files_labelled_ABCDEFGHI) we will now filter out trees with poor topologies according to a boostrap test (n=1000) performed with [IQ-TREE](http://www.iqtree.org). 

In order to do that we first need to reformat the input. We'll do these operations in folder [08_bootstrapping_check](./08_bootstrapping_check), which contains some sample files. One of the steps involved a script from https://github.com/vinuesa/get_phylomarkers

We'll start by creating pruned FASTA files with all diploids and only one allopolyploid sequence for each aligned cluster:

```
# Remove sequences with only gaps using trimAl

for FILE in *corrected_simplified.sort.fna; do
echo $FILE;
trimal/source/trimal -in $FILE -out $FILE.trimmed.fna -noallgaps;
done

# Extract multifasta using list of references for each species and label (alleles)
# Example for all labels (A, B, C, D, E, F, G, H and I) of B422 species. Repeat in each files and species
# See file Bmex_A.list.txt

for FILE in *.trimmed.fna; do cat Bmex_A.list.txt | \
	awk '{gsub("_","\\_",$0);$0="(?s)^>"$0".*?(?=\\n(\\z|>))"}1' | pcregrep -oM -f - $FILE > $FILE.extract_Bmex_A.fna; done

# Find and save FASTA files including allopoliploid sequences
# Example: codes of allopolyploids species: Bmex; Bret; Bboi; Bhyb; Brup; Bpho and B422

grep -E 'Bmex|Bret|Bboi|Bhyb|Brup|Bpho|B422' *.fna | cut -d":" -f1 > list_files_extract_allopoly.txt


# copy and save those files:

ls | grep -f list_files_extract_allopoly.txt | xargs mv -t DIRECTORY


# Create non-paremetric bootstrapping trees using IQ-TREE
# Input corresponds to FASTA alinged files with diploid sequences plus polyploid sequences.
# Example: 99716_c37819_extract_Bmex_A.fna
# Output: 1000 bootstrapping trees for each fna file

ls *.fna | parallel --gnu -j 50 iqtree-1.6.9-Linux/bin/iqtree -b 1000 -nt 1 -AICc -s {} :::
```

Now we are ready to re-label the bootstrapped trees in folder [iqtree_non_parametric_bootstrap_1000](././08_bootstrapping_check/iqtree_non_parametric_bootstrap_1000), producing relabelled files:
```
perl scripts/_bootstrap_label_stats.pl iqtree_non_parametric_bootstrap_1000

# Example output: 98131_c38723_extract_Bmex_A.fna.boottrees.relabelled
```

Now, as earlier, we should check the diploid backbone of all bootstrapped, relabelled and rooted trees:
```
for FILE in *.root.ph; do
perl scripts/_check_diploids.pl $FILE Osat;
done &>> stats_diploid_skeleton_1000_boot.log

# print summaries of bootstrap topologies, with bootstrapped trees numbered from 1 to 1000
awk 'c&&!--c;/# .\//{print $0; c=8}' stats_diploid_skeleton_1000_boot.log | \
	cut -d" " -f1,2 | sed 's/# .\///g' | awk '{printf "%s%s",$0,NR%2?"\t":RS}' | \
	sort -t/ -k1,1 -k2,2n > table_files_topologies_1000.tsv

# Example:
# 100546_c28312.extract_B422_E.fna.boottrees.relabelled/1.root.ph    Osat,Hvul,Bsta,Bdis,Barb,Bpin,Bsyl,


# Find acepted topologies, imposed by user
# 'Osat,Hvul,Bsta,Bdis,Barb,Bsyl,Bpin,'; 'Osat,Hvul,Bsta,Bdis,Barb,Bpin,Bsyl,'

sed -r 's/\/[0-9]+\.root.ph//g' table_files_topologies_1000.tsv | grep -e 'Osat,Hvul,Bsta,Bdis,Barb,Bsyl,Bpin,' -e 'Osat,Hvul,Bsta,Bdis,Barb,Bpin,Bsyl,' | uniq -c > STATS_acepted_topologies_1000_bootstrapping.txt


# sum up topologies for the same files:

awk '{a[$2]+=$1}END{for(i in a) print i,a[i]}' STATS_acepted_topologies_1000_bootstrapping.txt

# Keep names in each files with correct topology

cat table_files_topologies_1000.tsv | grep -e 'Osat,Hvul,Bsta,Bdis,Barb,Bsyl,Bpin,' -e 'Osat,Hvul,Bsta,Bdis,Barb,Bpin,Bsyl,' 

# New trimAl job to remove columns with all gaps after filtering

for FILE in *.fna; do
echo $FILE;
trimal/source/trimal -in $FILE -out $FILE.trimmed.fna -noallgaps -keepseqs;
done

# Concatenate genes partitions

ls *.trimmed.fna > list.txt
get_phylomarkers/concat_alignments.pl list.txt > MSA.fna
```

The next step in this section would be to remove underrepresented labels, which in our case we chose to be those present in less than 12 (6%) gene trees. 


## 6) Consensus subgenome labelled trees and MSA files

In this step the allopolyploid allele labels that we had been using (lower case) were simplified to consensus subgenomes labels (in capitals, as in the paper). This was guided by a four-way analysis, described in the paper, which included:

+ Patristic distances and their clustering in a Principal Component Analysis (see Excel file [patristic_distances.xlsx](./09_consensus_labels/patristic_distances.xlsx) produced by [Geneious](https://www.geneious.com/) from the concatenated ML tree obtained in the previous step. 

+ Allelle percentages in bootstrap trees (see Excel file [bootstrap_label_percentages.xlsx](./09_consensus_labels/bootstrap_label_percentages.xlsx)). 

+ Relative frequency of allele labels per species (up to four alleles).

+ The topology of the ML tree. 

The final subgenomic labelling rules are:

```
allele   consensus sungenome label
a+c             A
b               B
d               D
e               E
e(Bret)         E_core
f+g+h+i         H

```

Finally, sequences for each species and allele were first collapsed and then consensus computed with https://github.com/josephhughes/Sequence-manipulation/blob/master/Consensus.pl:
```
perl -lne 'if(/^(>.*)/){ $head=$1 } else { $fa{$head} .= $_ } END{ foreach $s (keys(%fa)){ print "$s\n$fa{$s}\n" if($s =~ /Bpho_F/ || $s =~ /Bpho_G/ || $s =~ /Bpho_H/|| $s =~ /Bpho_I/) }}' MSA.fasta > Bpho_F_G_H_I.fasta

perl Sequence-manipulation/Consensus.pl -iupac -in Bpho_F_G_H_I.fasta -out Bpho_H_consensus.fasta
```

## Data files of figures in the paper

The files containing sequence alignments, trees and XML config files for [BEAST](https://beast.community) and cross-bracing analyses will be available after publication at [https://doi.org/10.5061/dryad.ncjsxksqw](https://doi.org/10.5061/dryad.ncjsxksqw)

























