#!/usr/bin/env Rscript 

# This script implements the subgenome assignment algorithm
#
# Example call: ./homeolog_distribution.R 10 INPUT_homeolog_distribution.tsv output.txt

# Tested in R (version 3.5.1; 2018-07-02)
# By: Antonio Díaz-Pérez (june 03, 2020)

# FORMAT OF TEXT INPUT FILE:
#***************************

# The first row must indicate the column names

# For row 2,...,(n+1), where "n" is the number of genes:

#column 1 = ID_allele
#column 2 = label given to allele according to the "Nearest Diploid Species Node" algorithm
#column 3,4,...,(j+2) = number of bootstrap replicates grafted in each of "j" 
#                       branches of the diploid species tree

args = commandArgs(trailingOnly=TRUE)

if(length(args) < 3) {
	stop("# usage: ./homeolog_distribution.R <perc cutoff> <infile> <outfile>.n", call.=FALSE)
} else {
	threshold.1 = as.numeric( args[1] ) # threshold value in % for cluster construction
	input.name = args[2] 
	output.name = args[3]
}


threshold <- threshold.1/100
allele.matrix <- read.table(input.name,header=T)
number.branches <- length(allele.matrix[1,])-2
allele.factor <- as.factor(allele.matrix[,2])
number.alleles <- length(levels(allele.factor))
name.alleles <- levels(allele.factor)
branch.names <- colnames(allele.matrix)[3:(2+number.branches)]

for(i in name.alleles) {

	sub.matrix <- allele.matrix[which(allele.matrix[,2] == i),]
	branch.vector <- vector()
	for(j in 1:number.branches) {
		branch.vector[j] <- sum(sub.matrix[,(2+j)])
	}

	max.value <- branch.vector[which.max(branch.vector)]
	higher.threshold <- which(branch.vector >= max.value*threshold)
	higher.names <- branch.names[higher.threshold]
	higher.values <- branch.vector[higher.threshold] 
	higher.values.2 <- as.matrix((higher.values/max.value)*100)
	higher.values.3 <- higher.values.2

	for(k in 1:length(higher.values.2[,1])) {
		higher.values.3[k,1] <- paste(format( higher.values.2[k,1],
			nsmall=2,digits=2),collapse=" ") 
   }

	# append results to output file 
	write(" ",file=output.name,append=T)
	write(paste("Bootstrap frequency per branch for allele ",i,":",sep=""),
		file=output.name,append=T)
	write.table(cbind(branch.names,branch.vector),file=output.name,append=T,
		quote = F, row.names = F, col.names = F)

	write(" ",file=output.name,append=T)
	write(paste("Branches with bootstrap support equal/higher than ",
		threshold*100,"% (threshold) of the maximum value:",sep=""),
		file=output.name,append=T)
	write.table(cbind(as.matrix(higher.names),as.matrix(higher.values.3)),
		file=output.name,append=T,quote = F, row.names = F, col.names = F)
   write(" ",file=output.name,append=T)
   write("---------------------------------------------------------",
		file=output.name,append=T)
}

