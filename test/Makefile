# Makefile for testing complete pipeline

cwd=${PWD}
scripts_dir="${cwd}/../scripts"

outdir01="${cwd}/01_collapsed_core_clusters_rename_toy"
outdir02="${cwd}/02_blocks_Brachypodium_toy"

# Filter alignments to get diploid blocks
# input: ${outdir01} output: ${outdir02}
nearest_diploid_species_node_1:
	mkdir -p ${outdir02} ;
	for FILE in ${outdir01}/*.fa ; do \
        BSFILE="$(basename $$FILE)" ; \
		echo ${BSFILE} ; \
		#perl ${scripts_dir}/_trim_MSA_block.pl -i $$FILE -o ${outdir02}/${BSFILE}.block.fna -m 100 -M 100 -O 0.5 -R 'c\d+_g\d+_i' -t ${scripts_dir}/species2rename.tsv &>> ${outdir01}/blocks_Brachypodium_toy.log ; \
    done



#nearest_diploid_species_node_1:
#	for FILE in ${outdir01}/*.fa ; do \ 
#		echo ${FILE} ; \ 
#	done


#mkdir -p ${outdir02};
#perl ${scripts_dir}/_trim_MSA_block.pl -i ${FILE} -o ${outdir02}/${FILE}.block.fna \ 
#			-m 100 -M 100 -O 0.5 -R 'c\d+_g\d+_i' ${scripts_dir}/species2rename.tsv &>> \ 
#			${outdir01}/blocks_Brachypodium_toy.log; \ 
#	done 
#	@echo Save list of files (block.fna) with all taxa
#	grep -e '^>.*_Osat' -e '^>.*_Hvul' -e '^>.*_Bsta' -e '^>.*_Bdis' \
#		-e '^>.*_Barb' -e '^>.*_Bpin' -e '^>.*_Bsyl' ${outdir02}/*.block.fna | cut -d":" -f1 | \
#		uniq -c | cut -d" " -f7,8 | grep "^7" | cut -d" " -f2 > ${outdir02}/list_block_all_fasta_toy.txt
