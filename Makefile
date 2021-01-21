# Makefile for installation and test

cwd=${PWD}

## third-party dependencies

# trimal, see https://github.com/scapella/trimal
trimdir="trimal-1.4.1"
trimtar="v1.4.1.tar.gz"
trimalurl="https://github.com/scapella/trimal/archive/${trimtar}"

# newick_utils, see http://cegg.unige.ch/newick_utils
nudir="newick-utils-1.6"
nutar="newick-utils-1.6-Linux-x86_64-enabled-extra.tar.gz"
nuurl="http://cegg.unige.ch/pub/${nutar}"

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

install:
	make trimal newick_utils iqtree get_homologues concat_aln consensus

test:
	perl allopolyploids.t

trimal:
	@echo "Download and compile ${trimurl}"
	cd ${cwd}/bin; wget -c ${trimalurl}; tar xfz ${trimtar}; cd ${trimdir}/source; make

newick_utils:
	@echo "Download and compile ${nuurl}"
	cd ${cwd}/bin; wget -c ${nuurl}; tar xfz ${nutar}; cd ${nudir}; ./configure; make

iqtree:
	@echo "Download ${iqurl}"
	cd ${cwd}/bin; wget -c ${iqurl}; tar xfz ${iqtar}

get_homologues:
	@echo "Download and install ${ghurl}"
	cd ${cwd}/bin; git clone ${ghurl}; cd ${ghdir}; perl install.pl force

concat_aln:
	@echo "Download ${caurl}"
	cd ${cwd}/bin; wget -c ${caurl}

consensus:
	@echo "Download ${courl}"
	cd ${cwd}/bin; wget -c ${courl}

