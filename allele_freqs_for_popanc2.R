##This is a script that formats the output of entropy for the program popanc
#this script was written for A. Springer C. maculatus admixture data
##this file needs to have the following format: 


##########################################################################################

# The popanc software requires three input text files that contain the genetic data: one file
# for each of the source (parental) populations and one for the admixed population.
# Genetic data files: Each genetic data file provides the genotypic data for one population.
# A header row in each file gives the number of SNPs and the number of individuals separated
# by a space. This is followed by one row per SNP that begins with the chromosome (or
#                                                                                  scaffold) number and the position (in centi-Morgans or other units) separated by a colon.
# After the locus information, the genotypic data is given for each individual as the count
# of the non-reference allele. This need not be an integer value, if, for example, genotypes
# were imputed or estimated from low or moderate coverage DNA sequence data. Here is an
# example with 5 SNPs and 10 individuals.
#
# 5 10
# 1:20 0 0 0 1 1 0 0 0 0 1
# 1:30 2 2 1 1 2 1 0 0 1 2
# 1:35 0 0 0 1 0 1 0 0 0 1
# 1:42 0 1 0 2 1 1 2 0 0 0
# 1:55 1.8 1.2 0 2 1.2 0.9 0.1 1.9 2 2

##########################################################################################


#the input file for this script is a list of scaffold names and positions
#for every SNP in the clean, filtered SNP set
#in this case, I have 79,078 SNPs in my set

#this input file was created by using GREP to pull out just the scaffold names 
#and positions from the .gl file of genotype likelihoods
#the scaffold name and position are the first item in each line of the .gl file

#first, read in file of scaffold names/positions for each SNP:
joe <- read.csv("AS_cmac_admix_ng_SNP_positions.txt", stringsAsFactors = FALSE, header = FALSE)

#check that nrow matches the number of SNPs in ourvariant set: 

# nrow(joe)
# [1] 79079

##########################################################################################

##next, read in the genotype matrix: 

#number of loci (rows): 
L <- 79079
#number of individuals (columns):
N <- 1536

G <- matrix(scan("AS_cmac_admix_MCMC_summaries_gprob_k3_NOHEADERS", n = N*L, sep = ","), nrow = N, ncol = L, byrow = T)

#Need to transpose G in order to be correct for popanc... need individuals as columns, not rows
G <- t(G)

##########################################################################################

#now get the IDs for individuals so we can sort the matrix by population etc.

##read in IDs file: 
#TO AVOID ANY ISSUES WITH THINGS BEING SORTED DIFFERENTLY IN THE .GL FILE, DRAW THE IDS DIRECTLY FROM THE .GL FILE: 
ids <- as.matrix(read.csv("cmac_admix_gl_order_ids.txt", header = FALSE, sep = " ", stringsAsFactors = FALSE))
#make into vector:
ids <- ids[1,]

#fix the one mislabeled BZ rep:

ids[grep(ids, pattern = "BZ.5.F5.L")] <- sub(pattern = "BZ.5.F5.L", 
                                                ids[grep(ids, pattern = "BZ.5.F5.L")], 
                                                replacement = "BZ.5.F20.L")

#and fix the four individuals labeled F20_21:

ids[grep(ids, pattern = "F20_21.L")] <- sub(pattern = "F20_21.L", 
                                            ids[grep(ids, pattern = "F20_21.L")], 
                                            replacement = "F20.L")



##########################################################################################

#read in my 3 allele frequency files for BF, BZ, and CA on cowpea (C) for the F1 generation: 
setwd("./estpEM_allele_freqs_newgenome")

BF_freq <- read.csv("AS.BF.1.F1.C_estpEM_allele_freqs.txt", sep = " ", stringsAsFactors = FALSE, header = FALSE)
BZ_freq <- read.csv("AS.BZ.1.F1.C_estpEM_allele_freqs.txt", sep = " ", stringsAsFactors = FALSE, header = FALSE)
CA_freq <- read.csv("AS.CA.1.F1.C_estpEM_allele_freqs.txt", sep = " ", stringsAsFactors = FALSE, header = FALSE)

setwd("../")

##now separate out just the SNPs that best inform ancestry for each population pair:

#absolute difference in allele frequencies between parental pops desired: 
abs_diff <- 0.2

#Here are the columns that need to be kept for each treatment group:
SNPs_for_BF_BZ <- which(abs(BF_freq[,3] - BZ_freq[,3]) > abs_diff)
SNPs_for_BF_CA <- which(abs(BF_freq[,3] - CA_freq[,3]) > abs_diff)
SNPs_for_BZ_CA <- which(abs(BZ_freq[,3] - CA_freq[,3]) > abs_diff)


##########################################################################################

#now need to separate out just the beetles from each treatment group
#and sort out the columns required for each treatment group

#the grep below sorts out beetles by pop, rep, generation, and host
#for some reason dots don't work well in R grep, but oh well

#grep(ids, pattern = "BF.[1].F[12][890].L")
#also need to add the number of SNPs (nrow) and number of individuals (ncol)

#first, parents: 
ids[grep(ids, pattern = "BF.[1].F1.C")]


BF_for_BFxBZ_popanc <- cbind(joe[SNPs_for_BF_BZ,],
                             G[SNPs_for_BF_BZ, grep(ids, pattern = "BF.[1].F1.C")])
BF_for_BFxBZ_popanc <- rbind(c(nrow(BF_for_BFxBZ_popanc), 
                               ncol(BF_for_BFxBZ_popanc), 
                               rep("", times = ncol(BF_for_BFxBZ_popanc) - 2 )), 
                               BF_for_BFxBZ_popanc)
write.table(BF_for_BFxBZ_popanc, file = "BF_for_BFxBZ_for_popanc", quote = FALSE,
                               row.names = FALSE, col.names = FALSE)


BZ_for_BFxBZ_popanc <- cbind(joe[SNPs_for_BF_BZ,],
                             G[SNPs_for_BF_BZ, grep(ids, pattern = "BZ.[1].F1.C")])
BZ_for_BFxBZ_popanc <- rbind(c(nrow(BZ_for_BFxBZ_popanc), 
                               ncol(BZ_for_BFxBZ_popanc), 
                               rep("", times = ncol(BZ_for_BFxBZ_popanc) - 2 )), 
                             BZ_for_BFxBZ_popanc)
write.table(BZ_for_BFxBZ_popanc, file = "BZ_for_BFxBZ_for_popanc", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

BF_for_BFxCA_popanc <- cbind(joe[SNPs_for_BF_CA,],
                             G[SNPs_for_BF_CA, grep(ids, pattern = "BF.[1].F1.C")])
BF_for_BFxCA_popanc <- rbind(c(nrow(BF_for_BFxCA_popanc), 
                               ncol(BF_for_BFxCA_popanc), 
                               rep("", times = ncol(BF_for_BFxCA_popanc) - 2 )), 
                             BF_for_BFxCA_popanc)
write.table(BF_for_BFxCA_popanc, file = "BF_for_BFxCA_for_popanc", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

CA_for_BFxCA_popanc <- cbind(joe[SNPs_for_BF_CA,],
                             G[SNPs_for_BF_CA, grep(ids, pattern = "CA.[1].F1.C")])
CA_for_BFxCA_popanc <- rbind(c(nrow(CA_for_BFxCA_popanc), 
                               ncol(CA_for_BFxCA_popanc), 
                               rep("", times = ncol(CA_for_BFxCA_popanc) - 2 )), 
                             CA_for_BFxCA_popanc)
write.table(CA_for_BFxCA_popanc, file = "CA_for_BFxCA_for_popanc", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

BZ_for_BZxCA_popanc <- cbind(joe[SNPs_for_BZ_CA,],
                             G[SNPs_for_BZ_CA, grep(ids, pattern = "BZ.[1].F1.C")])
BZ_for_BZxCA_popanc <- rbind(c(nrow(BZ_for_BZxCA_popanc), 
                               ncol(BZ_for_BZxCA_popanc), 
                               rep("", times = ncol(BZ_for_BZxCA_popanc) - 2 )), 
                             BZ_for_BZxCA_popanc)
write.table(BZ_for_BZxCA_popanc, file = "BZ_for_BZxCA_for_popanc", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

CA_for_BZxCA_popanc <- cbind(joe[SNPs_for_BZ_CA,],
                             G[SNPs_for_BZ_CA, grep(ids, pattern = "CA.[1].F1.C")])
CA_for_BZxCA_popanc <- rbind(c(nrow(CA_for_BZxCA_popanc), 
                               ncol(CA_for_BZxCA_popanc), 
                               rep("", times = ncol(CA_for_BZxCA_popanc) - 2 )), 
                             CA_for_BZxCA_popanc)
write.table(CA_for_BZxCA_popanc, file = "CA_for_BZxCA_for_popanc", quote = FALSE,
            row.names = FALSE, col.names = FALSE)

##########################################################################################

#Finally, need to get the F7 and F20 comparisons for all admixed lineages
#there should be:
#BF/BZ 1:5, F7 and 20, C and L = 5 x 2 x 2 = 20
#BF/CA 1:5, F7 and 20, C and L = 5 x 2 x 2 = 20
#BZ/CA 1:5, F7 and 20, C and L = 5 x 2 x 2 = 20

#need to split out the ids to a more usable form: 
#recall that "." is a special character so needs to be escaped
treatments <- character(length = length(ids))

split_ids <- strsplit(ids, split = "\\.")
for(i in 1:length(ids)){
  treatments[i] <- paste(split_ids[[i]][1], 
                         split_ids[[i]][2], 
                         split_ids[[i]][3], 
                         split_ids[[i]][4], sep = ".")
}

#check that I have the right number of admixed populations: 
admixed_treatments <- unique(treatments[grep(treatments, pattern = "_")])
#should be 60 total

#make function that can be used in a loop for each pop type

make_popanc_input_file <- function(gl_file, SNPs_list, scaff_names, pattern){
  popanc_matrix <- cbind(scaff_names[SNPs_list,],
                         gl_file[SNPs_list, grep(ids, pattern = pattern)])
  popanc_matrix <- rbind(c(nrow(popanc_matrix), 
                       ncol(popanc_matrix), 
                       rep("", times = ncol(popanc_matrix) - 2 )), 
                       popanc_matrix)
  write.table(popanc_matrix, file = paste(pattern, "_for_popanc", sep = ""), 
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  print(paste(pattern, " nSNPs is ", nrow(popanc_matrix), 
                                          " & nInds is ", ncol(popanc_matrix)))
}

#make files for BF_BZ:
for(i in 1:20){
  make_popanc_input_file(gl_file = G, SNPs_list = SNPs_for_BF_BZ, 
                         scaff_names = joe, pattern = admixed_treatments[i])
}

#make files for BF_CA:
for(i in 21:40){
  make_popanc_input_file(gl_file = G, SNPs_list = SNPs_for_BF_CA, 
                         scaff_names = joe, pattern = admixed_treatments[i])
}

#make files for BZ_CA:
for(i in 41:60){
  make_popanc_input_file(gl_file = G, SNPs_list = SNPs_for_BZ_CA, 
                         scaff_names = joe, pattern = admixed_treatments[i])
}

