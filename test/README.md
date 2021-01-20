This folder is used to test the complete pipeline with a toy Brachypodium dataset.

The starting dataset is in folder 01_collapsed_core_clusters_rename_toy 
which contains 42 core sequence clusters that include all diploids and two 
polyploids (B.hybridum [Bhyb] and B.phoenicoides [Bpho]) taxa.

These were obtained as follows:

```
# 1.1) Run get_homologues-est
# see //eead-csic-compbio.github.io/get_homologues/manual-est/

nohup perl get_homologues/get_homologues-est.pl -d genome_transcripts \
	-m cluster -I genome_transcripts/no_sorghum_no_sylCor.list -M -A -S 80

# this produces
# genome_transcripts_est_homologues/arb8075_alltaxa_no_sorghum_no_sylCor.list_algOMCL_e0_S80_/
#
# from which a subset of 42 clusters was saved in folder 00_1_get_homologues_toy/ (not in repo)

# 1.2) Compute core clusters and save them in folder 00_2_core_clusters_Brachypodium_toy (not in repo)

perl get_homologues/compare_clusters.pl -d 00_1_get_homologues_toy \
	-o 00_2_core_clusters_Brachypodium_toy -m -n 

# 1.3) Collapse and save locally align clusters in 00_2_core_clusters_Brachypodium_toy (not in repo)

for FILE in `ls 00_2_core_clusters_Brachypodium_toy/*.fna`; do
	echo $FILE;
	perl get_homologues/annotate_cluster.pl -D -f $FILE -o $FILE.aln.fna -c 20
done

# 1.4) Simplify FASTA headers and shorten cluster file names

# Original header:
#>c31981_g1_i1_chr7:19409754:19415236:+[arb_80_75] collapsed:1,{6}, Pfam:Golgi complex component 7 (COG7);(PF10191,)
# simplified header:
#>c31981_g1_i1_chr7[arb_80_75]

for FILE in *collapsed.fna.aln.fna; do
	echo $FILE;
	perl -p -i -e 's/>(.+?) .+/>$1/g; s/:\d+:\d+:[+-]//g' $FILE;
done

for f in *collapsed.fna.aln.fna;
	do mv -- "$f" "$(basename -- "$f" .collapsed.fna.aln.fna).fa";
done

# The resulting files are those contained in folder 01_collapsed_core_clusters_rename_toy/
```
