# PhyloSD: Phylogenomic Subgenome Detection pipeline

[![Build Status](https://travis-ci.com/eead-csic-compbio/allopolyploids.svg?branch=master)](https://travis-ci.com/eead-csic-compbio/allopolyploids)

This pipeline was designed by Ruben Sancho, Antonio Diaz, Pilar Catalan and Bruno Contreras Moreira for the subgenome identification of homeologous diploid genomes present in allopolyploids, also considering potentially unknown progenitors. 

The Phylogenomic Subgenomic Detection (PhyloSD) pipeline includes three sequential Nearest Diploid Species Node, Bootstrapping Refinement, and Subgenome Assignment algorithms.

The protocol is explained in detail in the article (publication in process) titled:

PhyloSD: Phylogenomic detection of known and ghost subgenomes of polyploid plants
R Sancho, LA Inda, A Díaz-Pérez, DL Des Marais, SP Gordon, J Vogel, B Contreras-Moreira, P Catalán

A Docker container with all dependencies preinstalled can be found at https://hub.docker.com/repository/docker/eeadcsiccompbio/allopolyploids

<!-- made with perl -lne 'if(/^(#{1,}) (.*)/){ ($i,$t)=($1,$2); $l=lc($t); $l=~s/\W/\-/g; print "$i [$t](#$l)"}'-->

- [0) Pre-processing of the data set and installation of dependencies](#0-pre-processing-of-the-data-set-and-installation-of-dependencies)
	- [Data set and abbreviations](#data-set-and-abbreviations)
	- [Run GET_HOMOLOGUES-EST to construct the core gene/transcript alignments](#run-get_homologues-est-to-construct-the-core-genetranscript-alignments)
- [1) NEAREST DIPLOID SPECIES NODE algorithm](#1-nearest-diploid-species-node-algorithm)
- [2) BOOTSTRAPPING REFINEMENT algorithm](#2-bootstrapping-refinement-algorithm)
- [3) SUBGENOME ASSIGNMENT algorithm](#3-subgenome-assignment-algorithm)
- [4) OPTIONAL DATING ANALYSIS](#4-optional-dating-analysis)



## 0) Pre-processing of the data set and installation of dependencies

This pipeline has several dependencies, which can be installed as explained below in folder [bin/](./bin):

|software|URL|
|-------|---|
|TrimAl|https://github.com/scapella/trimal|
|newick_utils|http://cegg.unige.ch/newick_utils|
|IQ-TREE|http://www.iqtree.org|
|GET_HOMOLOGUES-EST|https://github.com/eead-csic-compbio/get_homologues|
|concat_alignments.pl|https://github.com/vinuesa/get_phylomarkers|
|Consensus.pl|https://github.com/josephhughes/Sequence-manipulation|
|R|https://www.r-project.org|

In addition, the instructions below require **wget**, **curl**, **make**, **git**, **perl**  and **parallel** on your system.

In Debian-like systems like Ubuntu please copy and paste the following command lines in your terminal to install:

```
# Perl
sudo apt-get install -y libdb-dev parallel git curl cpanminus
curl -L https://cpanmin.us | perl - --sudo App::cpanminus
cpanm --sudo -v --installdeps --notest --cpanfile cpanfile .

# R 
sudo apt-get install -y r-base

# third-party binaries
make install
```

In RedHat-like systems:
```
# Perl
sudo yum install libdb-devel
curl -L https://cpanmin.us | perl - --sudo App::cpanminus
cpanm --sudo -v --installdeps --notest --cpanfile cpanfile .

# R
sudo yum install epel-release
sudo apt-get install R

# third-party binaries
make install
```

In order to check you installation please run:
```
make test
```


### Data set and abbreviations

In order to download the data used for the *Brachypodium*, used in the tutorial below,
[benchmark](https://github.com/eead-csic-compbio/allopolyploids/releases/download/1.0/Brachypodium_bench.tar.gz) 
please do:
```
make brachy 
```

Note that *Brachypodium* RNA-seq data were also deposited at the 
[European Nucleotide Archive](https://www.ebi.ac.uk/ena) with run accessions numbers:

* ERR3633153 (B. arbuscula)
* ERR3634426 (B. boissieri)
* ERR3634464 (B. hybridum)
* ERR3634593 (B. mexicanum) *
* ERR3634717 (B. phoenicoides Bpho422)
* ERR3634894 (B. phoenicoides Bpho6)
* ERR3634970 (B. pinnatum)
* ERR3636695 (B. retusum)
* ERR3636791 (B. rupestre) 
* ERR3636828 (B. stacei)

*B. sylvaticum* RNA-seq data from accession Brasy-Esp were obtained from the study by 
[Fox et al., 2013](https://doi.org/10.3732/apps.1200011). 

Transcriptomic data were also retrieved for the outgroup species Oryza sativa (SRX738077) 
and Hordeum vulgare (ERR159679), which were used to root the phylogenetic *Brachypodium* trees 
and to provide a stem branch for the grafting of ancestral homeologs/subgenomes.

The following abbreviations were used:

```
Bmex --> B.mexicanum
Bret --> B.retusum
Bboi --> B.boissieri
Bhyb --> B.hybridum
Brup --> B.rupestre
Bpho --> B.phoenicoides (accession Bpho6)
B422 --> B.phoenicoides (accession Bpho422)
Bsta --> B.stacei
Bdis --> B.distachyon
Barb --> B.arbuscula
Bpin --> B.pinnatum
Bsyl --> B.sylvaticum
Osat --> Oryza sativa
Hvul --> Hordeum vulgare
```

### Run GET_HOMOLOGUES-EST to construct the core gene/transcript alignments

This first step requires https://github.com/eead-csic-compbio/get_homologues and the set of transcripts in folder genome_transcripts, most of them assembled de novo with https://github.com/trinityrnaseq/trinityrnaseq . With the following command we request clusters with sequence identity >= 80% and 75% coverage of the shortest sequence

```
bin/get_homologues/get_homologues-est.pl -d genome_transcripts \
	-m cluster -I genome_transcripts/species.list -M -A -S 80 \
	&> gen.M.A.S80.core.cluster.log
```
This produces 3675 clusters.

We can now compute pan-gene matrices for these core clusters:
```
bin/get_homologues/compare_clusters.pl \
-d genome_transcripts_est_homologues/arb8075_alltaxa_no_sorghum_no_sylCor.list_algOMCL_e0_S80_/ \
-o core_clusters_Brachypodium -m -n &> compare.core.log
```

We can optionally plot an Average Nucleotide Identity (ANI) matrix:

```
# Note: this requires additional R dependencies
bin/get_homologues/install_R_deps.R
bin/get_homologues/plot_matrix_heatmap.sh -i \
   genome_transcripts_est_homologues/arb8075_alltaxa_no_sorghum_no_sylCor.list_algOMCL_e0_S80_Avg_identity.tab \
   -H 14 -W 26 -t "ANI of transcripts in 3675 core clusters" -N -o pdf
```

From the clusters obtained earlier we can now produce MVIEW collapsed, multiple sequence alignments (MSA), while also annotating Pfam domains:
```
for FILE in `ls core_clusters/*.fna`; do
   echo $FILE;
	bin/get_homologues/annotate_cluster.pl -D -f $FILE -o $FILE.aln.fna -c 20 &>> collapse.align.core.log;
done
```


## 1) NEAREST DIPLOID SPECIES NODE algorithm


### 1.1) Edit polyconfig.pm file to adapt the previous information to your specific analyses. Define diploids, polyploids, outgroup, rooted species, ad-hoc labelling rules, ...


### 1.2) Simplify headers and names of files from collapsed core clusters

We might need to simplify the FASTA headers, which is something that should be taylored according to the user's data.

```
for FILE in *collapsed.fna.aln.fna; do
echo $FILE;
perl -p -i -e 's/>(.+?) .+/>$1/g; s/:\d+:\d+:[+-]//g' $FILE;
done

Original header:
>c43592_g1_i1_chr8:11969473:11975162:+[B422_80_75] collapsed:6,{12}, Pfam:Microtubule associated protein (MAP65/ASE1 family);(PF03999,)

Simplified header:
>c43592_g1_i1_chr8[B422_80_75]
```
The goal is to leave the species identifier as the last bit, in this case within [brackets]. You might need to follow a sligthly different approach with your own data.

You can also simplify the names of your files (optional):
```
for f in *.fna;
do mv -- "$f" "$(basename -- "$f" .collapsed.fna.aln.fna).fa";
done

Original names:
99998_c33211_g1_i1_chr2-26015115-26017843--.collapsed.fna.aln.fna

Renamed names:
99998_c33211_g1_i1_chr2-26015115-26017843--.fa
```

**Note:** All FASTA files are compressed in fasta.files.tar.bz2 and saved in 01_collapsed_core_clusters directory


### 1.3) Filter Mutiple Sequences Alignments (MSAs) according to the diploid block and remove inconsistent positions

We now take these alignments of nucleotide sequences of both diploid and polyploid species and produce trimmed FASTA files suitable for phylogenetic tree inference. The goal is to define a solid diploid backbone, which should be covered by outgroup sequences as well, and then use it to filter out polyploid sequences/alleles with diploid block overlap < $MINBLOCKOVERLAP. Therefore, some sequences and species might be removed from the initial input.

This process is computed using scripts/_trim_MSA_block.pl and the defauld parameters are defined in polyconfig.pm file:

Default values for paremeters controlling how aligned blocks of sequences are produced in script _trim_MSA_block.pl
our $MINBLOCKLENGTH = 100;
our $MAXGAPSPERBLOCK = 100; # tolerated gaps for diploids in block
our $MINBLOCKOVERLAP = 0.50; # fraction of diploid block covered by outgroups & polyploids

We run the tool:
```
for FILE in *.fa;
do echo $FILE;
perl scripts/_trim_MSA_block.pl -i $FILE -o $FILE.block.fna \
-m 100 -M 100 -O 0.5 -R 'c\d+_g\d+_i' \
-t scripts/species2rename.tsv &>> blocks.log;
done
```
We recover 2001 aligned core clusters with a consistent diploid block. You can save the filter files in 02_blocks directory

**Warning**: In this step, some alignments can lost outgroups (Osat or Hvul). If we count the number of block.fna files with each taxon we recover:
Osat --> 1925/2001
Hvul --> 1944/2001
Bsta --> 2001/2001
Bdis --> 2001/2001
Barb --> 2001/2001
Bpin --> 2001/2001
Bsyl --> 2001/2001

You can count the files with all diploid and outgroup taxa using the regular expression:
```
grep -e '^>.*_Osat' -e '^>.*_Hvul' -e '^>.*_Bsta' -e '^>.*_Bdis' -e '^>.*_Barb' -e '^>.*_Bpin' -e '^>.*_Bsyl' *.block.fna \
| cut -d":" -f1 | uniq -c | cut -d" " -f7,8 | grep -c "^7"
```
We recover 1878 alignments with all diploid and outgroup taxa.

We save a list of files with all diploid and ougroup taxa:
```
grep -e '^>.*_Osat' -e '^>.*_Hvul' -e '^>.*_Bsta' -e '^>.*_Bdis' -e '^>.*_Barb' -e '^>.*_Bpin' -e '^>.*_Bsyl' *.block.fna | cut -d":" -f1 | uniq -c | cut -d" " -f7,8 | grep "^7" | cut -d" " -f2 > list_block_all_fasta.txt
```

We can now trim the resulting blocks with https://vicfero.github.io/trimal to remove the inconsistent positions of the alignments:
```
for FILE in *block.fna; do
echo $FILE;
bin/trimal/source/trimal -in $FILE -out $FILE.trimmed.fna -automated1;
done
```


### 1.4) Compute ML (Maximun-Likelihood) gene trees

We can now run IQ-TREE in parallel. The results are saved in the folder 03_iqtree_treefiles/iqtree_treefile (see results in treefile.tar.bz2 file)
```
ls *.trimmed.fna | parallel --gnu -j 3 bin/iqtree-1.6.12-Linux/bin/iqtree \
   -alrt 1000 -bb 1000 -nt AUTO -AICc -s {} :::
```


### 1.5) Root and sort nodes in trees, loosing bootstrap and aLRT supports on the way

We will now root and ladderize the nodes (starting with rice. Root defined in polyconfig.pm) in the resulting trees, which are stored in folder 03_iqtree_treefiles/rooted_trees (see root.ph.tar.bz2 file).

**Note:** this script requires Perl modules Bio::TreeIO and Bio::Phylo:
```
sudo cpan -i Bio::TreeIO Bio::Phylo

for FILE in *.treefile;
do perl scripts/_reroot_tree.pl $FILE > $FILE.root.ph;
echo $FILE;
done
```

**Warning:** if there is not "root" taxon, an empty file is created. Remove this files are recomended
```
find . -empty -type f -delete
```
**Warning:** some files have not Hvul outgroup. This will be filtered upstream. We save them into a list.

We save the list of files (root.ph) with all taxa:
```
ls *.ph | grep -f ../02_blocks/list_block_all_fasta.txt > list_root_tree_all_taxa.txt
```
**Note:** one of this files would be lost in the next step because of "Calling node's branch length will be zero (set -FORCE to force)--exiting". This file starts with "12615".


### 1.6) Check diploid skeleton (topology) for each tree

We can now ask how many different diploid skeletons are there and select gene trees with congruent diploid topology. The next script must be run within folder 03_iqtree_treefiles/rooted_trees and requires http://cegg.unige.ch/newick_utils
```
for FILE in *root.ph;
do perl scripts/_check_diploids.pl $FILE;
done > diploids.log
```
The tool also create the pruned tree files (trees with only diploid species. Not save in this example).

**Note:** There are 1924 pruned trees (and 1925 rooted trees) because one tree show the error "Calling node's branch length will be zero (set -FORCE to force)--exiting". This file starts with "12615"

We check and count the diferent diploid topologies:
```
grep -v "#" diploids.log | sort  | uniq -c | sort -n
```

Examples (number of trees/genes and diploid topologies):
e.g.
Counts Diploid skeleton/topologies
 75    (Osat,(Hvul,(Bsta,(Bdis,(Barb,(Bsyl,Bpin))))));
124    (Osat,(Hvul,(Bsta,(Bdis,(Bsyl,(Barb,Bpin))))));
165    (Osat,(Hvul,(Bsta,(Bdis,(Bpin,(Barb,Bsyl))))));
254    (Osat,(Hvul,(Bsta,(Bdis,(Barb,(Bpin,Bsyl))))));

Congruent diploid skeletons/topologies (all variants):

(Osat,(Hvul,(Bsta,(Bdis,(Barb,(Bpin,Bsyl))))));
(Osat,(Hvul,(Bsta,(Bdis,(Barb,(Bsyl,Bpin))))));

(Osat,(Hvul,(Bsta,(Bdis,((Bsyl,Bpin),Barb)))));
(Osat,(Hvul,(Bsta,(Bdis,((Bpin,Bsyl),Barb)))));

We convert the .log file to a table (tsv format):
```
grep -e '# .*fna.treefile.root.ph*' -e '(' "diploids.log" \
| cut -d" " -f1,2 | sed 's/# //g' | awk '{printf "%s%s",$0,NR%2?"\t":RS}' > table_diploids.tsv
```

We create the table with all congruent diploid topologies:
```
grep -w \
-e '(Osat,(Hvul,(Bsta,(Bdis,(Barb,(Bpin,Bsyl))))));' \
-e '(Osat,(Hvul,(Bsta,(Bdis,(Barb,(Bsyl,Bpin))))));' \
-e '(Osat,(Hvul,(Bsta,(Bdis,((Bsyl,Bpin),Barb)))));' \
-e '(Osat,(Hvul,(Bsta,(Bdis,((Bpin,Bsyl),Barb)))));' \
table_diploids.tsv  > table_diploids_congruent.tsv
```
We check and count the congruent diploid topologies
```
cut -f2 table_diploids_congruent.tsv | sort | uniq -c
```
Counts Diploid skeleton/topologies
254    (Osat,(Hvul,(Bsta,(Bdis,(Barb,(Bpin,Bsyl))))));
 75    (Osat,(Hvul,(Bsta,(Bdis,(Barb,(Bsyl,Bpin))))));

We recover 329 trees with the congruent diploid topology.

We create the list of congruent FASTA files (save in folder 02_blocks).
```
cut table_diploids_congruent.tsv -f1 | sed 's/.treefile.root.ph//' > list_congruent_fasta.txt
```
We create the list of congruent tree (rooted.ph) files (save in folder 03_iqtree_treefiles/rooted_trees).
```
cut table_diploids_congruent.tsv -f1 > list_congruent_trees.txt
```


### 1.7) We will now copy files (.fna, root.ph, root.ph.pruned) with congruent diploid topologies to folder 04_congruent_and_labelled_files

Run the commands in 02_blocks
```
ls | grep -f list_congruent_fasta.txt | xargs cp -t ../04_congruent_and_labelled_files
```
Run in 03_iqtree_treefiles/rooted_trees
```
ls | grep -f list_congruent_trees.txt | xargs cp -t ../../04_congruent_and_labelled_files
```
**Note**: These results are not included but you can run the code to recover them


### 1.8) Labelling polyploid homeologs

In this step, we will label polyploid sequences (homeolog) with some pre-defined codes according to their position with respect to the diploid skeleton. We will do that with scripts/_check_lineages_polyploids.pl, which has several parameters defined previously (see polyconfig.pm file).
This script takes as input an MSA and tree file with congruent diploid skeletons and produces two labelled MSA files with sorted diploid species and labelled polyploid homeologs (A,B,C,etc) species, a full MSA with some only-gap rows where one polypoid sequence is missed and a reduced MSA excluding missing polyploid sequences (homeologs), and the labelled trees. Ad-hoc rules encoded (see polyconfig.pm) are defined using the species or clades of diploids used such as a reference to do the labelling process. These rules should be carefully adapted to other clades. Labelled files are stored in the folder 04_congruent_and_labelled_files with the names label.fna.tar.bz2 (full labelled MSAs), label.reduced.fna.tar.bz2 (reduced labelled MSAs) and label.ph.tar.bz2 (Labelled trees), respectively.

```
for FILE in *.fna; do
perl scripts/_check_lineages_polyploids.pl -t $FILE\.treefile.root.ph -f $FILE;
done
```
You can run the script using a loop:
```
nohup loop_check_lineage_polyploid.sh &> lineage_codes_blocks_ABCDEFGHI.log &
```

The preliminary statistic of labelled polyploid homeologs can be calculated using the scripts/_make_lineages_stats.pl.
**Note:** these labels correspond to homeologs, which are represented as lower-case letters in the phylogenies of the paper.
```
perl scripts/_make_lineages_stats.pl lineage_codes_blocks_ABCDEFGHI_Brachypodium.log
```
```
Table of preliminary homeolog-type statistics:

        A       B       C       D       E       F       G       H       I       all
Bmex    111     54      51      7       16      0       1       0       0       240
Bboi    96      48      58      12      37      3       7       1       1       263
Bret    46      13      49      6       83      22      40      18      15      292
Bhyb    2       157     2       114     0       0       0       0       1       276
Brup    1       1       1       5       62      39      84      48      57      298
Bpho    3       0       5       3       80      53      72      34      50      300
B422    1       2       2       4       70      53      79      34      51      296
     																									-------
   																				   Total alleles: 1965
```



## 2) BOOTSTRAPPING REFINEMENT algorithm

We use 05_bootstrapping_analysis such as working directory of this step


### 2.1) Set the pruned FASTA alignments (diploids + outgroups + one polyploid homeolog)

We create the pruned aligment using such as input files the reduced labelled MSAs (04_congruent_and_labelled_files/label.reduced.fna.tar.bz2).
Firstly, we do the files that will be include in the pruned alignment. One list/file per polyploid homeolog (see folder 05_bootstrapping_analysis/list_taxa)

For each labelled homeolog (above e.g. Bmex_A):
```
for FILE in *.label.reduced.fna.trimmed.fna; do cat ./list_taxa_Brachypodium/Bmex_A.list.txt \
| awk '{gsub("_","\\_",$0);$0="(?s)^>"$0".*?(?=\\n(\\z|>))"}1' | pcregrep -oM -f - $FILE > $FILE.extract_Bmex_A.fna; done
```

Some pruned MSAs (.extract_Bmex_A.fna, ...) will not include any polyploid. Thus, we have to save only the pruned MSA with polyploid sequences.
```
grep -E 'Bmex|Bboi|Bret|Bhyb|Brup|Bpho|B422' *.fna | cut -d":" -f1 > list_files_extract_with_polyploid.txt
```

Move pruned FASTA files with diploids + outgroups + one homeolog to one_allopolyploid_plus_diploids/2_1_files_extract_with_polyploid directory (now empty)
```
ls *.fna | grep -f list_files_extract_with_polyploid.txt | xargs mv -t files_extract_with_polyploid
```
Removed the other pruned FASTA ".extract" files without polyploid sequences


### 2.2) Run 1000 non-parametric bootstrapping replicates (iqtree)

```
ls *.fna | parallel --gnu -j 60 bin/iqtree-1.6.12-Linux/bin/iqtree -b 1000 -nt 1 -AICc -s {} :::
```
You can run it using a bash script and nohup:
```
nohup iqtree_non_parametric_bootstrapping_1000.sh > nohup_bootstrapping_1000.txt &
```

These .boottrees files will be saved in the folder 2_2_iqtree_boottrees_1000.


### 2.3) Extract, root and sort the bootstrap trees (boottrees) from the 1000 bootstrapping files

The script/_bootstrap_label_stats.pl is used for this task. If only root (rootonly) flag is actived, two files for each boottree are created: trees (.ph) extracted from .boottrees, and rooted and sorted trees (.root.ph)
```
nohup scripts/_bootstrap_label_stats.pl 2_2_iqtree_boottrees_1000/ rootonly &
```
The statistics have to take in count just in trees with the congruent diploid skeleton (see next step).


### 2.4) Check diploid skeletons of extracted, rooted and sorted boottrees (1000 per allele)

```
for FILE in ./*.relabelled/*.root.ph;
do perl scripts/_check_diploids.pl $FILE;
done &>> stats_diploid_skeleton_1000boottrees.log
```

See diploid skeletons topologies in "stats_diploid_skeleton_1000boottrees.log" created previously 

grep -e '# .*fa.block.label.reduced.fna.trimmed*' -e '(' "stats_diploid_skeleton_1000boottrees.log" \
| cut -d" " -f1,2 | sed 's/# .\///g' | awk '{printf "%s%s",$0,NR%2?"\t":RS}' \
| sort -t/ -k1,1 -k2,2n > table_diploid_skeleton_1000boottrees.tsv

1,965,000 boottrees (1,965 alleles/homeolog x 1,000 boottrees per homeolog),

Create a file with the list of directories (alleles/homeologs):
```
grep "# ./" stats_diploid_skeleton_1000boottrees.log | cut -d"/" -f2 \
| uniq > stats_diploid_skeleton_1000boottrees_directories.txt
```

Again, we check the congruent diploid topology (all variants) and save the files with this diploid skeleton:

(Osat,(Hvul,(Bsta,(Bdis,(Barb,(Bpin,Bsyl))))));
(Osat,(Hvul,(Bsta,(Bdis,(Barb,(Bsyl,Bpin))))));

(Osat,(Hvul,(Bsta,(Bdis,((Bsyl,Bpin),Barb)))));
(Osat,(Hvul,(Bsta,(Bdis,((Bpin,Bsyl),Barb)))));

We create a table with the files and all congruent topologies:
```
grep -w \
-e '(Osat,(Hvul,(Bsta,(Bdis,(Barb,(Bpin,Bsyl))))));' \
-e '(Osat,(Hvul,(Bsta,(Bdis,(Barb,(Bsyl,Bpin))))));' \
-e '(Osat,(Hvul,(Bsta,(Bdis,((Bsyl,Bpin),Barb)))));' \
-e '(Osat,(Hvul,(Bsta,(Bdis,((Bpin,Bsyl),Barb)))));' \
table_diploid_skeleton_1000boottrees.tsv  > table_diploid_skeleton_1000boottrees_congruent.tsv
```
A total of 1,002,542 boottrees show the congruent diploid skeleton

We create a table with the congruent directories:
```
cut -d"/" -f1 table_diploid_skeleton_1000boottrees_congruent.tsv \
| uniq > table_diploid_skeleton_1000boottrees_congruent_directories.txt
```
1,960 directories/alleles/homeologs show the congruent diploid skeleton

**Note:* five directories do not appear in this list because there are not any boottrees with congruent diploid skeleton. However, those should be considered to remove them as part of the "congruent_diploid_skeleton_1000" files with less than 100 boottrees showing congruent diploid skeleton.
```
diff stats_diploid_skeleton_1000boottrees_directories.txt table_diploid_skeleton_1000boottrees_congruent_directories.txt
```

Five directories (alleles) have 0 boottrees with congruent diploid skeleton:
```
108136_c31693_g1_i1_chr5-23779430-23783492-+.fa.block.label.reduced.fna.trimmed.fna.extract_Bboi_B.fna.boottrees.relabelled
108136_c31693_g1_i1_chr5-23779430-23783492-+.fa.block.label.reduced.fna.trimmed.fna.extract_Bpho_E.fna.boottrees.relabelled
111530_c38754_g2_i5_chr5-9236742-9244175-+.fa.block.label.reduced.fna.trimmed.fna.extract_Bhyb_B.fna.boottrees.relabelled
93237_c39035_g2_i1_chr2-30560610-30570642-+.fa.block.label.reduced.fna.trimmed.fna.extract_Bret_E.fna.boottrees.relabelled
93349_c39377_g2_i1_chr4-4428242-4434064--.fa.block.label.reduced.fna.trimmed.fna.extract_Bboi_C.fna.boottrees.relabelled
```

If you execute the folowing command, it returns nothing because there are not any cogruent diploid skeletons:
```
grep -A8 "111530_c38754_g2_i5_chr5-9236742-9244175-+.fa.block.label.reduced.fna.trimmed.fna.extract_Bhyb_B.fna.boottrees.relabelled" \
stats_diploid_skeleton_1000boottrees.log \
| grep -w -e '(Osat,(Hvul,(Bsta,(Bdis,(Barb,(Bpin,Bsyl))))));' \
-e '(Osat,(Hvul,(Bsta,(Bdis,(Barb,(Bsyl,Bpin))))));' \
-e '(Osat,(Hvul,(Bsta,(Bdis,((Bsyl,Bpin),Barb)))));' \
-e '(Osat,(Hvul,(Bsta,(Bdis,((Bpin,Bsyl),Barb)))));'
```

We check and count the congruent diploid topologies (1000 boottrees)
```
cut -f2 table_diploid_skeleton_1000boottrees_congruent.tsv | sort | uniq -c
```
```
Counts Diploid skeletons/topologies
929020 (Osat,(Hvul,(Bsta,(Bdis,(Barb,(Bpin,Bsyl))))));
 73522 (Osat,(Hvul,(Bsta,(Bdis,(Barb,(Bsyl,Bpin))))));
```

So, there are 1002542 boottrees of total 1965000 with congruent diploid skeleton.
We count the boottrees per each allele/homeolog with congruent diploid skeleton.
```
sed -r 's/\/[0-9]+\.root.ph//g' table_diploid_skeleton_1000boottrees_congruent.tsv | uniq -c \
| awk '{a[$2]+=$1}END{for(i in a) print i,a[i]}' > table_diploid_skeleton_1000boottrees_congruent_counts.tsv
```


### 2.5) Check files that show less than 100 bootrees with congruent diploid skeleton

We create a list of congruent files (boottrees):
```
cut -f1 table_diploid_skeleton_1000boottrees_congruent.tsv | sort -t"/" -k1,1 -k2,2n \
> list_diploid_skeleton_boottrees_congruent_files_sorted.txt
```
```
e.g. The list shows the directory/allele (100065_c35185..../) and the boottrees files with congruent diploid skeleton (1.root.ph, ...)

100065_c35185_g1_i1_chr2-2924216-2928564--.fa.block.label.reduced.fna.trimmed.fna.extract_Bboi_A.fna.boottrees.relabelled/1.root.ph
100065_c35185_g1_i1_chr2-2924216-2928564--.fa.block.label.reduced.fna.trimmed.fna.extract_Bboi_A.fna.boottrees.relabelled/2.root.ph
100065_c35185_g1_i1_chr2-2924216-2928564--.fa.block.label.reduced.fna.trimmed.fna.extract_Bboi_A.fna.boottrees.relabelled/3.root.ph
100065_c35185_g1_i1_chr2-2924216-2928564--.fa.block.label.reduced.fna.trimmed.fna.extract_Bboi_A.fna.boottrees.relabelled/4.root.ph
100065_c35185_g1_i1_chr2-2924216-2928564--.fa.block.label.reduced.fna.trimmed.fna.extract_Bboi_A.fna.boottrees.relabelled/5.root.ph
100065_c35185_g1_i1_chr2-2924216-2928564--.fa.block.label.reduced.fna.trimmed.fna.extract_Bboi_A.fna.boottrees.relabelled/6.root.ph
100065_c35185_g1_i1_chr2-2924216-2928564--.fa.block.label.reduced.fna.trimmed.fna.extract_Bboi_A.fna.boottrees.relabelled/7.root.ph
...
```

We have to check if any allele shows less than 100 boottrees with congruent diploid skeleton. In our analysis, 177 alleles (directories) show less than 100 boottrees with congruente diploid skeleton (five alleles show zero; see previous step).

e.g. Above you can see the directory/allele id and the counts of boottrees with congruent diploid skeleton
```
Allele id                                                                                                                    boottrees
108136_c31693_g1_i1_chr5-23779430-23783492-+.fa.block.label.reduced.fna.trimmed.fna.extract_Bhyb_D.fna.boottrees.relabelled   	1
110196_c13070_g1_i1_chr6-8389172-8390724-+.fa.block.label.reduced.fna.trimmed.fna.extract_Bpho_H.fna.boottrees.relabelled     	1
108136_c31693_g1_i1_chr5-23779430-23783492-+.fa.block.label.reduced.fna.trimmed.fna.extract_Bmex_C.fna.boottrees.relabelled   	2
110196_c13070_g1_i1_chr6-8389172-8390724-+.fa.block.label.reduced.fna.trimmed.fna.extract_B422_H.fna.boottrees.relabelled     	2
93349_c39377_g2_i1_chr4-4428242-4434064--.fa.block.label.reduced.fna.trimmed.fna.extract_Brup_F.fna.boottrees.relabelled      	2
110196_c13070_g1_i1_chr6-8389172-8390724-+.fa.block.label.reduced.fna.trimmed.fna.extract_Bhyb_D.fna.boottrees.relabelled     	3
...
```
We create a list (e.g. list_diploid_skeleton_boottrees_congruent_directories.txt) of congruent directories (alleles) with at least 100 boottrees with congruent diploid skeleton. This list includes 1,788 alleles (directories) (we do not include the 177 alleles (directories) with less than 100 boottrees with congruent diploid skeleton)
e.g. to sort the alleles of the list:
```
cat list_diploid_skeleton_boottrees_congruent_directories.txt | sort > list_diploid_skeleton_boottrees_congruent_directories_sorted.txt

100065_c35185_g1_i1_chr2-2924216-2928564--.fa.block.label.reduced.fna.trimmed.fna.extract_Bboi_A.fna.boottrees.relabelled
100065_c35185_g1_i1_chr2-2924216-2928564--.fa.block.label.reduced.fna.trimmed.fna.extract_Bhyb_D.fna.boottrees.relabelled
100065_c35185_g1_i1_chr2-2924216-2928564--.fa.block.label.reduced.fna.trimmed.fna.extract_Bmex_A.fna.boottrees.relabelled
100065_c35185_g1_i1_chr2-2924216-2928564--.fa.block.label.reduced.fna.trimmed.fna.extract_Bpho_E.fna.boottrees.relabelled
100128_c35206_g3_i1_chr2-21577376-21581870-+.fa.block.label.reduced.fna.trimmed.fna.extract_B422_F.fna.boottrees.relabelled
100128_c35206_g3_i1_chr2-21577376-21581870-+.fa.block.label.reduced.fna.trimmed.fna.extract_Bboi_E.fna.boottrees.relabelled
...
```

We create a list with the congruent files (boottrees) according to congruent directories (alleles) with at least 100 boottrees with congruent diploid skeleton:
```
grep -f list_diploid_skeleton_boottrees_congruent_directories_sorted.txt \
list_diploid_skeleton_boottrees_congruent_files_sorted.txt \
> list_diploid_skeleton_boottrees_congruent_files_filtered.txt
```


### 2.6) Copy congruent filtered files (boottrees) in 2_3_iqtree_boottrees_1000_congruent directory

We can run this command in the folder 2_2_iqtree_boottrees_1000
```
cat list_diploid_skeleton_boottrees_congruent_files_filtered.txt | \
xargs cp --parents -t ../2_3_iqtree_boottrees_1000_congruent/
```
You can check the number of boottrees into 2_3_iqtree_boottrees_1000_congruent directory
```
ls -R | grep -c "root.ph"

We count 993,743 boottrees with congruent diploid skeleton (1,788 directories, alleles)
```


### 2.7) Random selection of 100 boottrees (root.ph) for each directory/allele in the folder 2_3_iqtree_boottrees_1000_congruent

We run the commands in the folder 2_3_iqtree_boottrees_1000_congruent
```
for d in *.relabelled;
do for f in *.root.ph;
do shuf -n 100 -e $d/$f | xargs cp --parents -t ../2_4_iqtree_boottrees_100_random_congruent;
done; done
```

We can check the number of boottrees in 2_4_iqtree_boottrees_100_random_congruent directory´
```
ls -R | grep -c "root.ph"

We count 178,800 boottrees (1,788 directories/alleles x 100 boottrees per directory)
```

We count the number of boottrees per directory (allele) to make sure the number is 100:
```
for d in *.relabelled;
do (cd $d && ls | wc | cut -d" " -f5);
done

100
100
100
...
```

We put together in a single boottrees file each of the hundred trees of each allele (it will be create into each directory the boottrees file):
```
for d in *.relabelled;
do for f in *.ph;
do (cd $d && cat $f | sed 's/;/;\n/g' > $d.boottrees);
done; done
```

We move the new boottrees files to a new directory 2_5_iqtree_boottrees_100_random_congruent_concatenate (run in the folder 2_4_iqtree_boottrees_100_random_congruent)
```
for d in *.relabelled;
do for f in *.boottrees;
do (cd $d && mv $f ../../2_5_iqtree_boottrees_100_random_congruent_concatenate);
done; done
```

We can rename boottrees files in 2_5_iqtree_boottrees_100_random_congruent_concatenate directory:
```
for file in *.boottrees;
do mv "$file" "$(basename "$file" .boottrees.relabelled.boottrees).boottrees";
done
```


### 2.8) Re-label and compute the statistics of the 100 random boottrees with congruent diploid skeleton for each directory/allele:

We run the scripts/_bootstrap_label_stats.pl tool to re-label the boottrees and calculate the statistics of homeolog-types (A, B, C,...)
```
nohup scripts/_bootstrap_label_stats.pl \
2_5_iqtree_boottrees_100_random_congruent_concatenate/ > stats_100_boottrees_ramdom_congruent.tsv &
```


### 2.9) Additional round of re-labelling congruent boottrees

There may be alleles where not all of your boottrees have been labelled. This can occur if one polyploid homeolog branch is inserted between the two outgroups (g. Osat and Hvul) and you have not been defined any ad-hoc labelling rules for this possibility because the genus under study is more recent than both outgroups. To get more consistent statistics in this cases, we re-label all boottrees with a congruent diploid skeleton (not just 100 boottrees) and conduct the random selection of 100 re-labelled boottrees. 

Firstly, we create a list (e.g. list_alleles_round_2.txt in the folder 2_3_iqtree_boottrees_1000_congruent) of alleles (directories) that do not show 100 random boottrees re-labelled. We detect 79 alleles with some homeolog without labelling. We copy these to the filder 2_6_iqtree_boottrees_1000_congruent_round_2 and repeat the previous steps.
Run in the folder 2_3_iqtree_boottrees_1000_congruent
```
cat list_alleles_round_2.txt | xargs cp -r --parents -t ../2_6_iqtree_boottrees_1000_congruent_round_2/
```


### 2.10) We put together in a single boottrees file each of the 100 boottrees (additional round) of each allele.

Run in the folder 2_6_iqtree_boottrees_1000_congruent_round_2
```
for d in *.relabelled;
do for f in *.ph;
do (cd $d && cat $f | sed 's/;/;\n/g' > $d.boottrees);
done; done
```
We move the concatenated boottrees files to new directory 2_7_iqtree_boottrees_1000_congruent_round_2_concatenate
```
for d in *.relabelled;
do for f in *.boottrees;
do (cd $d && mv $f ../../2_7_iqtree_boottrees_1000_congruent_round_2_concatenate);
done; done
```
We rename the boottrees files in 2_7_iqtree_boottrees_1000_congruent_round_2_concatenate directory
```
for file in *.boottrees;
do mv "$file" "$(basename "$file" .boottrees.relabelled.boottrees).boottrees";
done
```


### 2.11) Re-label and compute the statistics of the boottrees (additional round) with congruent diploid skeleton for each directory/allele

We run the scripts/_bootstrap_label_stats.pl tool to re-label the boottrees and calculate the statistics of homeolog-types (A, B, C,...)
```
nohup scripts/_bootstrap_label_stats.pl \
2_7_iqtree_boottrees_1000_congruent_round_2_concatenate > stats_all_boottrees_congruent_round_2.tsv &
```


### 2.12) Random selection of 100 re-labelled boottrees (additional round) for each directory/allele

Run in 2_7_iqtree_boottrees_1000_congruent_round_2_concatenate to copy to 2_8_iqtree_boottrees_100_random_congruent_round_2 directory
```
for d in *.relabelled;
do for f in *.root.ph;
do shuf -n 100 -e $d/$f | xargs cp --parents -t ../2_8_iqtree_boottrees_100_random_congruent_round_2;
done; done
```
We check the amount of boottrees in 2_8_iqtree_boottrees_100_random_congruent_round_2 directory
```
ls -R | grep -c "root.ph"

7,900 boottrees (79 directories/alleles x 100 boottrees per directory)
```
We count the number of boottrees per directory (allele) to make sure the number is 100:
```
for d in *.relabelled;
do (cd $d && ls | wc | cut -d" " -f5);
done

100
100
100
...
```

we put together in a single boottrees file each of the hundred trees of each allele (it will be create into each directory the boottrees file)
```
for d in *.relabelled;
do for f in *.ph;
do (cd $d && cat $f | sed 's/;/;\n/g' > $d.boottrees);
done; done
```

We move the new boottrees files to the new 2_9_iqtree_boottrees_100_random_congruent_round_2_concatenate directory
```
for d in *.relabelled;
do for f in *.boottrees;
do (cd $d && mv $f ../../2_9_iqtree_boottrees_100_random_congruent_round_2_concatenate);
done; done
```

We recover the 79 .boottrees files and rename the boottrees files
```
for file in *.boottrees;
do mv "$file" "$(basename "$file" .boottrees.relabelled.boottrees).boottrees";
done
```


### 2.13) Re-label and compute the statistics of the random selection of 100 boottrees (additional round) with congruent diploid skeleton for each directory/allele. We run the scripts/_bootstrap_label_stats.pl tool to re-label the boottrees and calculate the statistics of homeolog-types (A, B, C,...)

```
nohup scripts/_bootstrap_label_stats.pl \
2_9_iqtree_boottrees_100_random_congruent_round_2_concatenate/ > stats_100_boottrees_ramdom_congruent_round_2.tsv &
```


### 2.14) Refinement of the alignments according to bootstrapping analysis

According to the congruent diploid skeleton statistics and the re-labelling of boottrees (two rounds), we replace with gaps those homeologs that do not reach the minimum threshold of 100 boottrees with right diploid skeleton, or wrong labelling according to the bootstrapping analyses. In some cases, the complete gene (with multiple homeolgos) with every homeolgo wrong can be removed. We copy the files "*.fa.block.label.fna.trimmed.fna (04_congruent_and_labelled_files)" to the nwe directory  "06_labelled_alignment_corrected". These files (congruent genes) are edited changing wrong labelled alleles, according to bootstrapping analyses, by gaps. Remember to update the labelled statistics subtracting the gapped alleles.
Finally 560 homeologs and 7 complete genes (alignment) are removed/gapped. We keep 1,405 homeologs and 322 alignments/genes (without removing under-represeted alleles (see next steps)). These file are saved such as compressed file fasta.label.trimmed.tar.bz2.


### 2.15) Run trimAl to remove columns with all gaps, kepping all gapped sequences

During the previous steps, some columns could have been transformed into "allgaps" columns. To check this, trimAl is used. We run the commands in 06_labelled_alignment_corrected. These files are saved such as compressed file fasta.label.trimmed.noallgaps.tar.bz2.

```
for FILE in *.fna; do
echo $FILE;
bin/trimal/source/trimal -in $FILE -out $FILE.noallgaps.fna -noallgaps -keepheader -keepseqs;
done
```


### 2.16) Remove the underrepresented alleles/homeolgos (keep only the homeologs with 10% of representation in the polyploid species) from .trimmed files (see /06_labelled_alignment_corrected/fasta.label.trimmed.tar.bz2)

We removed the following underrepresented alleles/homeologs for each gene (MSAs)
```
Bmex_D; Bmex_E; Bmex_F; Bmex_G; Bmex_H; Bmex_I
Bboi_D; Bboi_F; Bboi_G; Bboi_H; Bboi_I
Bret_B; Bret_D; Bret_F; Bret_H; Bret_I
Bhyb_A; Bhyb_C; Bhyb_E; Bhyb_F; Bhyb_G; Bhyb_H; Bhyb_I
Brup_A; Brup_B; Brup_C; Brup_D
Bpho_A; Bpho_B; Bpho_C; Bpho_D;
B422_A; B422_B; B422_C; B422_D;
```

We make a list of taxa includign dipoloids (Osat, Hvul, ...) + polyploid homeologs (Bmex_A, ..) represented (list_taxa.txt) in the new folder 07_labelled_alignment_corrected_filtered.
```
e.g.
Osat
Hvul
Bsta
Bdis
Barb
Bpin
Bsyl
Bmex_A
...
```
We copy the files .trimmed from 06_labelled_alignment_corrected directory to 07_labelled_alignment_corrected_filtered directory and execute the folowing lines:
```
for FILE in *.trimmed.fna; do cat list_taxa.txt \
| awk '{gsub("_","\\_",$0);$0="(?s)^>"$0".*?(?=\\n(\\z|>))"}1' | pcregrep -oM -f - $FILE > $FILE.filtered.fna; done
```
After this step, we check the MSA to confirm that there is still at least one homologous sequence.
```
grep -E 'Bmex|Bboi|Bret|Bhyb|Brup|Bpho|B422' *filtered.fna | cut -d":" -f1 | uniq > list_files_with_polyploid_allele_represented.txt
```
In this case we recover 322 MSA with at least one homeolgos, so all genes (MSAs)  have at least one polyploid allele (homeolog) and we do not remove any complete gen. These .filtered.fna files are saved such as compressed file fasta.label.trimmed.filtered.tar.bz2 in the folder 07_labelled_alignment_corrected_filtered


### 2.17) Check and remove "allgapped" columns

During the process some columns could have been transformed into "allgaps" columns. To check this, trimAl is used.

Run in 07_labelled_alignment_corrected_filtered directory

```
for FILE in *filtered.fna; do
echo $FILE;
bin/trimal/source/trimal -in $FILE -out $FILE.noallgaps.fna -noallgaps -keepheader -keepseqs;
done
```
The filtered and corrected files are saved such as compressed fasta.label.trimmed.filtered.noallgaps.tar.bz2 files.


### 2.18) Phylogenomic analysis of concatenated labelled, filtered  and corrected genes/MSAs ("Homeologs' ML consensus tree")

We create the coordinate file and concatenate the MSAs using get_phylomarkers/concat_alignments.pl tool
```
ls -1 *.filtered.fna.noallgaps.fna > list_genes.txt

bin/concat_alignments.pl list_genes.txt > MSA_labelled_filtered.fna
```
We format the concatenation_coordinates.txt to crate the partitional file input for iqtree execution.

```
grep -v -P "^#" concatenation_coordinates.txt | perl -lne 'if(/(\d+-\d+)$/){ print "DNA, part$. = $1" }' > MSA_labelled_filtered_partitions.txt
```
We run iqtree to infer the "Homeologs' ML consensus tree"
```
nohup bin/iqtree-1.6.12-Linux/bin/iqtree -s MSA_labelled_filtered.fna -spp MSA_labelled_filtered_partitions.txt -m MFP -AICc -alrt 1000 -bb 1000 -nt AUTO &
```

**Note:**This result (MSA_labelled_filtered_partitions.txt.treefile) is the "Homeologs' ML consensus tree" (in the paper the labelled homeologs are represented using lower-case letters (a, b,...)



## 3) SUBGENOME ASSIGNMENT algorithm

To find out if some of these homeologs (a, b, c ...) refer to the same subgenome, and therefore to amalgamate or collapse them, we use different approaches (chromosome counts, patristic distances, PCoA-MST (Principal coordinates analysis-minimum spanning tree), homeolog grafting distributions (from homeolog grafting frequencies from 100 re-labelled boottrees).
If the homeolog types match subgenomes (single-allelic subgenomes), it will be not necessary to amalgamate or collapse the homeolog sequences. However, if multiple homeolog types match one subgenome (compound-allelic subgenomes), it will be necessary to amalgamate these homeologs into one to obtain the corresponding subgenome.
According to chromosome counts, patristic distances, PCoA and homeolog distribution, we will compute the "Subgenomic ML consensus tree".
**Note:** The subgenomes are shown such as upper-case letters (A1, A2, B...) in the "Subgenomic ML consensus tree" of the paper.
**Note2:** The files of this stage are saved in the folder 08_labelled_alignment_corrected_filtered_subgenomes.


### 3.1) Chromosome counts

See previous studies and ourself chromosome counts and genomic size analyses (see paper).


### 3.2) Patristic distance

Pairwise patristic distances between the Brachypodium diploid orthologous branches and the polyploid homeologous subgenomic branches of the "Homeologs' ML consensus tree" based on 322 nuclear core genes, and 1,307 orthologous and homeologous sequences. Patristic distances were calculated with Geneious R11.1.5. See Excel file "08_1_patristic distances.xlsx"

**Note:**You can calculate the patristic distances using other softwares, includig open source softwares, such as R packages (e.g. adephylo and ape)


### 3.3) Principal Coordinate Analysis (PCoA) and superimposed Minimum Spanning Tree (MST) generated from patristic distances and the Homeologs’ ML consensus tree.
Separate Principal Coordinate Analysis (PCoA) was performed using the patristic distance values of the data matrices with NTSYS-pc v2.10j software, and minimum spanning trees (MST) were superimposed on each of the PCoA plots (Figures are shown in the paper).

**Note:**Alternative open source software can be used for this step.


### 3.4) Compute the homeolog grafting distributions

The homeolog grafting distributions are cumputed analysing the re-labelling bootstrapping frequencies (100 boottrees), removing underrepresented alleles, using the R-script script/homeolog_distribution.R.

```
Rscript scripts/homeolog_distribution.R 10 08_2_input_homeolog_distribution.tsv 08_2_output_homeolog_distribution.txt
```
**Note:** The value 10 corresponds to the 10% of the highest frequency value, which represents the lowest threshold allowed to include additional branches in the distribution of that homeolog-type


### 3.5) Amalgamate homeologs (a, b, ...) to infer the compound-allelic subgenomes (A1, A2, B, ...)
According to the previous analyses computed using the Subgenome Assignment algorithm, we can amalgamate or collapse the homeologs what correspond to the same subgenome. This step is developed in the folder 08_3_amalgamate_homeologs.
```
Correspondences between homeologs (lower-case) and subgenomes (upper-case):
Bhyb --> B (b); D (d)
Bmex --> A (a+b+c)
Bboi --> A (a+b+c+e)
Bret --> A (a+c); E_core (e+g)
Brup --> E (e); G (f+g+h+i)
Bpho, B422 --> E (e); G (f+g+h+i)
```

We extract from 07_labelled_alignment_corrected_filtered/MSA_labellel_filtered.fna the aligned sequences to collapse (consensus) for each compound-allelic subgenomes.

```
First, we extract the aligment including diploid, outgroup and single allelic subgenomes sequences:

perl -lne 'if(/^(>.*)/){ $head=$1 } else { $fa{$head} .= $_ } END{ foreach $s (keys(%fa)){ print "$s\n$fa{$s}\n" if($s =~ /Osat/ || $s =~ /Hvul/ || $s =~ /Bsta/ || $s =~ /Bdis/ || $s =~ /Barb/ || $s =~ /Bpin/ || $s =~ /Bsyl/ || $s =~ /Bhyb_B/ || $s =~ /Bhyb_D/ || $s =~ /Brup_E/ || $s =~ /Bpho_E/|| $s =~ /B422_E/) }}' MSA_labelled_filtered.fna > Diploids_outgroups_Bhyb_b_Bhyb_d_Brup_e_Bpho_e_B422_e.fna


Second, we extract the alignment for each compound-allelic subgenomes sequences:

perl -lne 'if(/^(>.*)/){ $head=$1 } else { $fa{$head} .= $_ } END{ foreach $s (keys(%fa)){ print "$s\n$fa{$s}\n" if($s =~ /Bmex_A/ || $s =~ /Bmex_B/ || $s =~ /Bmex_C/) }}' MSA_labelled_filtered.fna > Bmex_a_b_c.fna

perl -lne 'if(/^(>.*)/){ $head=$1 } else { $fa{$head} .= $_ } END{ foreach $s (keys(%fa)){ print "$s\n$fa{$s}\n" if($s =~ /Bboi_A/ || $s =~ /Bboi_B/ || $s =~ /Bboi_C/ || $s =~ /Bboi_E/) }}' MSA_labelled_filtered.fna > Bboi_a_b_c_e.fna

perl -lne 'if(/^(>.*)/){ $head=$1 } else { $fa{$head} .= $_ } END{ foreach $s (keys(%fa)){ print "$s\n$fa{$s}\n" if($s =~ /Bret_A/ || $s =~ /Bret_C/) }}' MSA_labelled_filtered.fna > Bret_a_c.fna

perl -lne 'if(/^(>.*)/){ $head=$1 } else { $fa{$head} .= $_ } END{ foreach $s (keys(%fa)){ print "$s\n$fa{$s}\n" if($s =~ /Bret_E/ || $s =~ /Bret_G/) }}' MSA_labelled_filtered.fna > Bret_e_g.fna

perl -lne 'if(/^(>.*)/){ $head=$1 } else { $fa{$head} .= $_ } END{ foreach $s (keys(%fa)){ print "$s\n$fa{$s}\n" if($s =~ /Bpho_F/ || $s =~ /Bpho_G/ || $s =~ /Bpho_H/|| $s =~ /Bpho_I/) }}' MSA_labelled_filtered.fna > Bpho_f_g_h_i.fna

perl -lne 'if(/^(>.*)/){ $head=$1 } else { $fa{$head} .= $_ } END{ foreach $s (keys(%fa)){ print "$s\n$fa{$s}\n" if($s =~ /B422_F/ || $s =~ /B422_G/ || $s =~ /B422_H/|| $s =~ /B422_I/) }}' MSA_labelled_filtered.fna > B422_f_g_h_i.fna

perl -lne 'if(/^(>.*)/){ $head=$1 } else { $fa{$head} .= $_ } END{ foreach $s (keys(%fa)){ print "$s\n$fa{$s}\n" if($s =~ /Brup_F/ || $s =~ /Brup_G/ || $s =~ /Brup_H/|| $s =~ /Brup_I/) }}' MSA_labelled_filtered.fna > Brup_f_g_h_i.fna
```

We can use the script Consensus.pl:
```
e.g.

perl bin/Consensus.pl -iupac -in Bmex_a_b_c.fna -out Bmex_A_consensus.fna
```

We concatenate the diploids, outgroups, single-allelic and compound-allelic (.consensus.fna files) subgenomes sequences. 
```
cat Diploids_outgroups_Bhyb_b_Bhyb_d_Brup_e_Bpho_e_B422_e.fna *consensus.fna > MSA_labelled_filtered_consensus.fna
```

We check sequences, rename and sort them, and we save the result such as MSA_labelled_filtered_consensus_final.fna in the folder 08_labelled_alignment_corrected_filtered_subgenomes.


### 3.6) Compute the Subgenomic ML consensus tree

**Note:** the coordinates of partitions are the same as previous step (07_labelled_alignment_corrected_filtered)

```
nohup bin/iqtree-1.6.12-Linux/bin/iqtree -s MSA_labelled_filtered_consensus_final.fna -spp MSA_labelled_filtered_consensus_partitions.txt -m MFP -AICc -alrt 1000 -bb 1000 -nt AUTO &
```
The tree is saved such as MSA_labelled_filtered_consensus_partitions.txt.treefile file.



## 4) OPTIONAL DATING ANALYSIS

Once we have inferred the subgenomic ML consensus tree, 
we can optionally perform dating analyses with [BEAST2](https://www.beast2.org), 
which should be installed accordingly (not provided). 
We have conducted this analysis in the folder 09_Beast2_analysis.

```
The final subgenomes inferred for each polyploid species are:

Bhyb --> B; D
Bmex --> A1
Bboi --> A2
Bret --> A2; E1
Brup --> E2; G
Bpho6, Bpho422 --> E2; G
```

We format the FASTA file MSA_labelled_filtered_consensus_final.fna to NEXUS file MSA_labelled_filtered_consensus_final.nex (saved such as compressed file nex.gz) and create the XML file MSA_labelled_filtered_consensus_final.xml (saved such as compressed file xml.gz) using BEAUTI software. Thus, we exectute Beast2:
```
DISPLAY="" nohup java -Xms200g -Xmx400g -jar beast/lib/beast.jar -beagle_SSE -instances 10 -threads 10 MSA_labelled_filtered_consensus_final.xml &
```

We check the log file to confirm that all parameters show ESS > 200.
```
java -Xms200g -Xmx300g -jar soft/Tracer_v1.6/lib/tracer.jar file.log
```

We run Treeannotator to compute the MCC tree.
```
nohup beast/bin/treeannotator -burnin 10 -heights median -lowMem MSA_labelled_filtered_consensus_final.trees MSA_labelled_filtered_consensus_beast2_phylogram.tre &
```
The dating MCC tree is shown in the file MSA_labelled_filtered_consensus_beast2_phylogram.tre
