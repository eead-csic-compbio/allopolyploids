# Makefile for testing complete pipeline

cwd=${PWD}
scripts_dir="${cwd}/../scripts"
bin_dir="${cwd}/../bin"

outdir01="${cwd}/01_collapsed_core_clusters_rename_toy"
outdir02="${cwd}/02_blocks_Brachypodium_toy"
outdir03="${cwd}/03_others_iqtree_files_Brachypodium_toy"
outdir04="${cwd}/04_iqtree_treefiles_rooted_Brachypodium_toy"

test:
	make nearest_diploid_species_node_1 nearest_diploid_species_node_2\
		nearest_diploid_species_node_3 nearest_diploid_species_node_4\
		nearest_diploid_species_node_5 clean

clean:
	rm -rf ${outdir02} ${outdir03} ${outdir04}

# Filter alignments to get diploid blocks
# input: ${outdir01} output: ${outdir02}
nearest_diploid_species_node_1:
	mkdir -p ${outdir02};
	for FILE in ${outdir01}/*.fa; do\
		BSFILE=$$(basename $$FILE);\
		perl ${scripts_dir}/_trim_MSA_block.pl -i $$FILE\
			-o ${outdir02}/$${BSFILE}.block.fna\
			-m 100 -M 100 -O 0.5 -R 'c\d+_g\d+_i'\
			-t ${scripts_dir}/species2rename.tsv\
			>> ${outdir01}/blocks_Brachypodium_toy.log;\
		echo ${outdir02}/$${BSFILE}.block.fna;\
	done;

# Remove inconsistent positions in the alignments (trimAl)
# input: ${outdir02} output: ${outdir02}
nearest_diploid_species_node_2:
	for FILE in ${outdir02}/*block.fna; do\
		BSFILE=$$(basename $$FILE);\
		${bin_dir}/trimal-1.4.1/source/trimal -in $$FILE\
			-out ${outdir02}/$${BSFILE}.trimmed.fna -automated1;\
	done

# Compute ML gene trees
# input: ${outdir02} output: ${outdir03} ${outdir04}
nearest_diploid_species_node_3:
	ls ${outdir02}/*.trimmed.fna | parallel --gnu -j 3\
		${bin_dir}/iqtree-1.6.12-Linux/bin/iqtree -nt AUTO -AICc -s {} ::: ;\
	mkdir -p ${outdir04};\
	mv ${outdir02}/*.treefile ${outdir04};\
	mkdir -p ${outdir03};\
    mv ${outdir02}/*.trimmed.fna.* ${outdir03};\

# Root trees with outgroup "Osat" (See polyconfig.pm)
# input: ${outdir04} output: ${outdir04}
nearest_diploid_species_node_4:
	for FILE in ${outdir04}/*.treefile; do\
		BSFILE=$$(basename $$FILE);\
		perl ${scripts_dir}/_reroot_tree.pl $$FILE > ${outdir04}/$${BSFILE}.root.ph;\
	done

# Check diploid skeleton (topology) for each tree/gen 
# input: ${outdir04} output: ${outdir04}
nearest_diploid_species_node_5:
	for FILE in ${outdir04}/*.root.ph; do\
		perl ${scripts_dir}/_check_diploids.pl $$FILE;\
	done > ${outdir04}/diploids_Brachypodium_toy.log;\
	tail ${outdir04}/diploids_Brachypodium_toy.log

