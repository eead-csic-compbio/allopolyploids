# Makefile for installation and test

cwd=${PWD}

.PHONY: test install

## third-party dependencies

# trimal, see https://github.com/inab/trimal
trimdir="trimal-1.4.1"
trimtar="v1.4.1.tar.gz"
trimalurl="https://github.com/inab/trimal/archive/refs/tags/${trimtar}"

# newick_utils, see https://github.com/tjunier/newick_utils
nudir="newick_utils"
nuurl="https://github.com/tjunier/${nudir}.git"

# IQ-TREE, see http://www.iqtree.org
iqtar="iqtree-1.6.12-Linux.tar.gz"
iqurl="https://github.com/Cibiv/IQ-TREE/releases/download/v1.6.12/${iqtar}"

# GET_HOMOLOGUES-EST, see https://github.com/eead-csic-compbio/get_homologues
ghdir="get_homologues"
ghurl="https://github.com/eead-csic-compbio/${ghdir}.git"

# concat_alignments, see https://github.com/vinuesa/get_phylomarkers
caurl="https://raw.githubusercontent.com/vinuesa/get_phylomarkers/master/concat_alignments.pl"

# consensus, see https://github.com/josephhughes/Sequence-manipulation
courl="https://raw.githubusercontent.com/josephhughes/Sequence-manipulation/master/Consensus.pl"

## Brachypodium benchmark data
brachytar="Brachypodium_bench.tar.gz"
brachyurl="https://github.com/eead-csic-compbio/allopolyploids/releases/download/1.0/${brachytar}"


install:
	make trimal newick_utils iqtree get_homologues concat_aln consensus

test:
	perl allopolyploids.t

trimal:
	@echo "Downloading and compiling ${trimurl}"
	cd ${cwd}/bin; wget -c ${trimalurl}; tar xfz ${trimtar}; cd ${trimdir}/source; make; cd ../..; rm ${trimtar}

newick_utils:
	@echo "Downloading ${nuurl}"
	cd ${cwd}/bin; git clone ${nuurl}; cd ${nudir}; autoreconf -fi; ./configure; make	

iqtree:
	@echo "Downloading ${iqurl}"
	cd ${cwd}/bin; wget -c ${iqurl}; tar xfz ${iqtar}; rm ${iqtar}

get_homologues:
	@echo "Downloading and installing ${ghurl}"
	cd ${cwd}/bin; git clone ${ghurl}; cd ${ghdir}; perl install.pl force

concat_aln:
	@echo "Downloading ${caurl}"
	cd ${cwd}/bin; wget -c ${caurl}

consensus:
	@echo "Downloading ${courl}"
	cd ${cwd}/bin; wget -c ${courl}

brachy:
	@echo "Downloading ${brachyurl}"
	wget -c ${brachyurl}; tar xfz ${brachytar}; rm ${brachytar}
