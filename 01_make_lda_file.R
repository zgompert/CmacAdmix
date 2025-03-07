## This is a script that takes a file of genotype point estimates (pntest_filename.txt file) 
## and produces a list of initial q values for use in running ENTROPY 

## read genotype estimates, i.e., g = 0 * Pr(g=0) + 1 * Pr(g=1) + 2 * Pr(g=2)
## row = locus, column = number of individuals
## replace L (locus) and N (number of individuals) below such that it matches your file: 

L <- 87050
N <- 1536

## read in your pntest_filename.txt file: 
G <- matrix(scan("pntest_morefilter_maxcov8000_mind0.filtered_minCov2x_maxND20_Chaydenii_variants.txt",
               n=N*L, sep=" "), nrow=L, ncol=N, byrow=T)

## pca on the genotype covariance matrix
pcgcov <- prcomp(x=t(G), center=TRUE, scale=FALSE)


## calculate kmeans and lda
library(MASS)

##calculate kmeans for k = 3 since I have 3 beetle parental populations:
k3 <- kmeans(pcgcov$x[,1:5], 3,iter.max=10, nstart=10, algorithm="Hartigan-Wong")

## calculate lda for k = 3:
ldak3 <- lda(x=pcgcov$x[,1:5], grouping=k3$cluster, CV=TRUE)


## write off the lda file for use in entropy: 
write.table(round(ldak3$posterior,5), file="cmac_admix_lda_k3.txt", quote=F, row.names=F ,col.names=F)


## when you run entropy use provide the input values as, e.g., -q cmac_admix_lda_k3.txt
## also set -s to something like 50
