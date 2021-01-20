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

install:
	make trimal newick_utils

test:
	perl allopolyploids.t

trimal:
	@echo "Download and compile ${trimurl}"
	cd ${cwd}/bin; wget -c ${trimalurl}; tar xfz ${trimtar}; cd ${trimdir}/source; make; cd ${cwd}

newick_utils:
	@echo "Download and compile ${nuurl}"
	cd ${cwd}/bin; wget -c ${nuurl}; tar xfz ${nutar}; cd ${nudir}; ./configure; make; cd ${cwd}
