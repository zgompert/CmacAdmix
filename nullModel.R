# null model, see also Rego et al. 2019
# uses PicMin to test for repeated evolution
library(PicMin)
library(scales)
library(xtable)

analyNullSimP<-function(p0=0.5,p1=0.5,t=100,ne=100,lower=TRUE){
    p0[p0<1/100]<-1/100
    p0[p0>99/100]<-99/100
    p1[p1<1/100]<-1/100
    p1[p1>99/100]<-99/100
    fst<-1-(1-1/(2*ne))^t
    alpha<-p0 * (1-fst)/fst
    beta<-(1-p0) * (1-fst)/fst
    pv<-pbeta(p1,alpha+0.001,beta+0.001,lower.tail=TRUE)
    pv<-apply(cbind(pv*2,(1-pv)*2),1,min)
    pv 
}

############### get data ##################

## list of allele frequency files
aff<-list.files(pattern="AS") ## 78 files
L<-79079 
snps<-read.table("ChSNPs.txt",header=FALSE)
head(snps)

## BF on lentil
fp0<-41
fpt<-42:45
K<-length(fpt)
aff[fpt]
#[1] "AS.BF.1.F19.L_estpEM_allele_freqs.txt"
#[2] "AS.BF.2.F20.L_estpEM_allele_freqs.txt"
#[3] "AS.BF.3.F19.L_estpEM_allele_freqs.txt"
#[4] "AS.BF.4.F20.L_estpEM_allele_freqs.txt"
ts<-c(19,19,19,19)
nes<-c(53.2,39.5,38.7,43.7)

p0<-read.table(aff[fp0],header=FALSE)
pvsBF<-matrix(NA,nrow=K,ncol=L)
for(k in 1:K){
	pt<-read.table(aff[fpt[k]],header=FALSE)
	pvsBF[k,]<-analyNullSimP(p0=p0[,3],p1=pt[,3],t=ts[k],ne=nes[k],lower=TRUE)
}

# get rid of MAF < 0.01 SNPs, too noisy, and non-LG
com<-which(p0[,3] > 0.01 & p0[,3] < 0.99 & is.na(snps[,1])==FALSE)
par(mfrow=c(2,2))
for(k in 1:K){
	plot(-log10(pvsBF[k,com]),pch=19)
}

## run PicMin
## generate correlation matrix for null, 10,000 is number of reps
n<-K
emp_p_null_dat <- t(replicate(10000, PicMin:::GenerateNullData(1.0, n, a=0.3, b=5, genes=length(com))))

# Calculate the order statistics' p-values for each simulation
emp_p_null_dat_unscaled <- t(apply(emp_p_null_dat ,1, PicMin:::orderStatsPValues))

# Use those p-values to construct the correlation matrix
null_pMax_cor_unscaled <- cor( emp_p_null_dat_unscaled )
null_pMax_cor_unscaled

# convert obs to ordered
opvsBF<-pvsBF[,com]
for(k in 1:K){
	opvsBF[k,]<-PicMin:::EmpiricalPs(opvsBF[k,],large_i_small_p=FALSE)
}

# run test and save results
numReps<-20000
picpBF<-rep(NA,length(com))
piceBF<-rep(NA,length(com))
for(i in 1:length(com)){
	cat(i,"\n")
	o<-PicMin:::PicMin(opvsBF[,i],null_pMax_cor_unscaled, numReps = numReps)
	picpBF[i]<-o$p
	piceBF[i]<-o$config_est
}

#summary(p.adjust(picpBF,method="fdr"))


## BZ on lentil
fp0<-67
fpt<-68:72
K<-length(fpt)
aff[fpt]
#[1] "AS.BZ.1.F19.L_estpEM_allele_freqs.txt"
#[2] "AS.BZ.2.F19.L_estpEM_allele_freqs.txt"
#[3] "AS.BZ.3.F19.L_estpEM_allele_freqs.txt"
#[4] "AS.BZ.4.F20.L_estpEM_allele_freqs.txt"
#[5] "AS.BZ.5.F20.L_estpEM_allele_freqs.txt"
ts<-c(19,19,19,19,19)
nes<-c(78.1,95.0,68.0,56.4,58.2)

p0<-read.table(aff[fp0],header=FALSE)
pvsBZ<-matrix(NA,nrow=K,ncol=L)
for(k in 1:K){
	pt<-read.table(aff[fpt[k]],header=FALSE)
	pvsBZ[k,]<-analyNullSimP(p0=p0[,3],p1=pt[,3],t=ts[k],ne=nes[k],lower=TRUE)
}

# get rid of MAF < 0.01 SNPs, too noisy, and non-LG
com<-which(p0[,3] > 0.01 & p0[,3] < 0.99 & is.na(snps[,1])==FALSE)
par(mfrow=c(3,2))
for(k in 1:K){
	plot(-log10(pvsBZ[k,com]),pch=19)
}

## run PicMin
## generate correlation matrix for null, 10,000 is number of reps
n<-K
emp_p_null_dat <- t(replicate(10000, PicMin:::GenerateNullData(1.0, n, a=0.3, b=5, genes=length(com))))

# Calculate the order statistics' p-values for each simulation
emp_p_null_dat_unscaled <- t(apply(emp_p_null_dat ,1, PicMin:::orderStatsPValues))

# Use those p-values to construct the correlation matrix
null_pMax_cor_unscaled <- cor( emp_p_null_dat_unscaled )
null_pMax_cor_unscaled

# convert obs to ordered
opvsBZ<-pvsBZ[,com]
for(k in 1:K){
	opvsBZ[k,]<-PicMin:::EmpiricalPs(opvsBZ[k,],large_i_small_p=FALSE)
}

# run test and save results
numReps<-20000
picpBZ<-rep(NA,length(com))
piceBZ<-rep(NA,length(com))
for(i in 1:length(com)){
	cat(i,"\n")
	o<-PicMin:::PicMin(opvsBZ[,i],null_pMax_cor_unscaled, numReps = numReps)
	picpBZ[i]<-o$p
	piceBZ[i]<-o$config_est
}

## CA on lentil
fp0<-73
fpt<-74:78
K<-length(fpt)
aff[fpt]
#[1] "AS.CA.1.F19.L_estpEM_allele_freqs.txt"
#[2] "AS.CA.2.F19.L_estpEM_allele_freqs.txt"
#[3] "AS.CA.3.F19.L_estpEM_allele_freqs.txt"
#[4] "AS.CA.4.F20.L_estpEM_allele_freqs.txt"
#[5] "AS.CA.5.F20.L_estpEM_allele_freqs.txt"
ts<-c(19,19,19,19,19)
nes<-c(109.3,87.7,97.5,85.2,97.0)

p0<-read.table(aff[fp0],header=FALSE)
pvsCA<-matrix(NA,nrow=K,ncol=L)
for(k in 1:K){
	pt<-read.table(aff[fpt[k]],header=FALSE)
	pvsCA[k,]<-analyNullSimP(p0=p0[,3],p1=pt[,3],t=ts[k],ne=nes[k],lower=TRUE)
}

# get rid of MAF < 0.01 SNPs, too noisy, and non-LG
com<-which(p0[,3] > 0.01 & p0[,3] < 0.99 & is.na(snps[,1])==FALSE)
par(mfrow=c(3,2))
for(k in 1:K){
	plot(-log10(pvsCA[k,com]),pch=19)
}

## run PicMin
## generate correlation matrix for null, 10,000 is number of reps
n<-K
emp_p_null_dat <- t(replicate(10000, PicMin:::GenerateNullData(1.0, n, a=0.3, b=5, genes=length(com))))

# Calculate the order statistics' p-values for each simulation
emp_p_null_dat_unscaled <- t(apply(emp_p_null_dat ,1, PicMin:::orderStatsPValues))

# Use those p-values to construct the correlation matrix
null_pMax_cor_unscaled <- cor( emp_p_null_dat_unscaled )
null_pMax_cor_unscaled

# convert obs to ordered
opvsCA<-pvsCA[,com]
for(k in 1:K){
	opvsCA[k,]<-PicMin:::EmpiricalPs(opvsCA[k,],large_i_small_p=FALSE)
}

# run test and save results
numReps<-20000
picpCA<-rep(NA,length(com))
piceCA<-rep(NA,length(com))
for(i in 1:length(com)){
	cat(i,"\n")
	o<-PicMin:::PicMin(opvsCA[,i],null_pMax_cor_unscaled, numReps = numReps)
	picpCA[i]<-o$p
	piceCA[i]<-o$config_est
}

## BF X BZ on C
fp0<-c(41,67)
fpt<-grep(pattern="AS.BF_BZ.[12345].F[0123456789]+.C",x=aff,fixed=FALSE)
K<-length(fpt)
aff[fpt]
# [1] "AS.BF_BZ.1.F20.C_estpEM_allele_freqs.txt"
# [2] "AS.BF_BZ.1.F7.C_estpEM_allele_freqs.txt" 
# [3] "AS.BF_BZ.2.F20.C_estpEM_allele_freqs.txt"
# [4] "AS.BF_BZ.2.F7.C_estpEM_allele_freqs.txt" 
# [5] "AS.BF_BZ.3.F20.C_estpEM_allele_freqs.txt"
# [6] "AS.BF_BZ.3.F7.C_estpEM_allele_freqs.txt" 
# [7] "AS.BF_BZ.4.F20.C_estpEM_allele_freqs.txt"
# [8] "AS.BF_BZ.4.F7.C_estpEM_allele_freqs.txt" 
# [9] "AS.BF_BZ.5.F20.C_estpEM_allele_freqs.txt"
#[10] "AS.BF_BZ.5.F7.C_estpEM_allele_freqs.txt" 
ts<-c(19,6,19,6,19,6,19,6,19,6)
nes<-c(177.9,177.9,180.1,180.1,169.7,169.7,163.3,163.3,128.0,128.0)

p0a<-read.table(aff[fp0[1]],header=FALSE)
p0b<-read.table(aff[fp0[2]],header=FALSE)
p0<-(p0a[,3]+p0b[,3])/2
pvsBFxBZ_C<-matrix(NA,nrow=K,ncol=L)
for(k in 1:K){
        pt<-read.table(aff[fpt[k]],header=FALSE)
        pvsBFxBZ_C[k,]<-analyNullSimP(p0=p0,p1=pt[,3],t=ts[k],ne=nes[k],lower=TRUE)
}

# get rid of MAF < 0.01 SNPs, too noisy, and non-LG
com<-which(p0 > 0.01 & p0 < 0.99 & is.na(snps[,1])==FALSE)
par(mfrow=c(3,2))
for(k in 1:6){
        plot(-log10(pvsBFxBZ_C[k,com]),pch=19)
}

## run PicMin
## generate correlation matrix for null, 10,000 is number of reps
n<-K/2
emp_p_null_dat <- t(replicate(10000, PicMin:::GenerateNullData(1.0, n, a=0.3, b=5, genes=length(com))))

# Calculate the order statistics' p-values for each simulation
emp_p_null_dat_unscaled <- t(apply(emp_p_null_dat ,1, PicMin:::orderStatsPValues))

# Use those p-values to construct the correlation matrix
null_pMax_cor_unscaled <- cor( emp_p_null_dat_unscaled )
null_pMax_cor_unscaled

# convert obs to ordered
opvsBFxBZ_C<-pvsBFxBZ_C[,com]
for(k in 1:K){
	opvsBFxBZ_C[k,]<-PicMin:::EmpiricalPs(opvsBFxBZ_C[k,],large_i_small_p=FALSE)
}

# run test and save results
numReps<-20000
picpBFxBZ_C_early<-rep(NA,length(com))
piceBFxBZ_C_early<-rep(NA,length(com))
picpBFxBZ_C_late<-rep(NA,length(com))
piceBFxBZ_C_late<-rep(NA,length(com))
early<-seq(2,10,2)
late<-seq(1,9,2)
for(i in 1:length(com)){
	cat(i,"\n")
	o<-PicMin:::PicMin(opvsBFxBZ_C[early,i],null_pMax_cor_unscaled, numReps = numReps)
	picpBFxBZ_C_early[i]<-o$p
	piceBFxBZ_C_early[i]<-o$config_est
	o<-PicMin:::PicMin(opvsBFxBZ_C[late,i],null_pMax_cor_unscaled, numReps = numReps)
	picpBFxBZ_C_late[i]<-o$p
	piceBFxBZ_C_late[i]<-o$config_est
}

## BF X BZ on L
fp0<-c(41,67)
fpt<-grep(pattern="AS.BF_BZ.[12345].F[0123456789]+.L",x=aff,fixed=FALSE)
K<-length(fpt)
aff[fpt]
# [1] "AS.BF_BZ.1.F18.L_estpEM_allele_freqs.txt"
# [2] "AS.BF_BZ.1.F5.L_estpEM_allele_freqs.txt" 
# [3] "AS.BF_BZ.2.F20.L_estpEM_allele_freqs.txt"
# [4] "AS.BF_BZ.2.F7.L_estpEM_allele_freqs.txt" 
# [5] "AS.BF_BZ.3.F20.L_estpEM_allele_freqs.txt"
# [6] "AS.BF_BZ.3.F7.L_estpEM_allele_freqs.txt" 
# [7] "AS.BF_BZ.4.F20.L_estpEM_allele_freqs.txt"
# [8] "AS.BF_BZ.4.F7.L_estpEM_allele_freqs.txt" 
# [9] "AS.BF_BZ.5.F20.L_estpEM_allele_freqs.txt"
#[10] "AS.BF_BZ.5.F7.L_estpEM_allele_freqs.txt" 
ts<-c(19,6,19,6,19,6,19,6,19,6)
nes<-c(85.7,85.7,54.1,54.1,64.4,64.4,72.9,72.9,74.1,74.1)

p0a<-read.table(aff[fp0[1]],header=FALSE)
p0b<-read.table(aff[fp0[2]],header=FALSE)
p0<-(p0a[,3]+p0b[,3])/2
pvsBFxBZ<-matrix(NA,nrow=K,ncol=L)
for(k in 1:K){
        pt<-read.table(aff[fpt[k]],header=FALSE)
        pvsBFxBZ[k,]<-analyNullSimP(p0=p0,p1=pt[,3],t=ts[k],ne=nes[k],lower=TRUE)
}

# get rid of MAF < 0.01 SNPs, too noisy, and non-LG
com<-which(p0 > 0.01 & p0 < 0.99 & is.na(snps[,1])==FALSE)
par(mfrow=c(3,2))
for(k in 1:6){
        plot(-log10(pvsBFxBZ[k,com]),pch=19)
}

## run PicMin
## generate correlation matrix for null, 10,000 is number of reps
n<-K/2
emp_p_null_dat <- t(replicate(10000, PicMin:::GenerateNullData(1.0, n, a=0.3, b=5, genes=length(com))))

# Calculate the order statistics' p-values for each simulation
emp_p_null_dat_unscaled <- t(apply(emp_p_null_dat ,1, PicMin:::orderStatsPValues))

# Use those p-values to construct the correlation matrix
null_pMax_cor_unscaled <- cor( emp_p_null_dat_unscaled )
null_pMax_cor_unscaled

# convert obs to ordered
opvsBFxBZ<-pvsBFxBZ[,com]
for(k in 1:K){
	opvsBFxBZ[k,]<-PicMin:::EmpiricalPs(opvsBFxBZ[k,],large_i_small_p=FALSE)
}

# run test and save results
numReps<-20000
picpBFxBZ_early<-rep(NA,length(com))
piceBFxBZ_early<-rep(NA,length(com))
picpBFxBZ_late<-rep(NA,length(com))
piceBFxBZ_late<-rep(NA,length(com))
early<-seq(2,10,2)
late<-seq(1,9,2)
for(i in 1:length(com)){
	cat(i,"\n")
	o<-PicMin:::PicMin(opvsBFxBZ[early,i],null_pMax_cor_unscaled, numReps = numReps)
	picpBFxBZ_early[i]<-o$p
	piceBFxBZ_early[i]<-o$config_est
	o<-PicMin:::PicMin(opvsBFxBZ[late,i],null_pMax_cor_unscaled, numReps = numReps)
	picpBFxBZ_late[i]<-o$p
	piceBFxBZ_late[i]<-o$config_est
}

## BF X CA on C
fp0<-c(41,73)
fpt<-grep(pattern="AS.BF_CA.[12345].F[0123456789]+.C",x=aff,fixed=FALSE)
K<-length(fpt)
aff[fpt]
# [1] "AS.BF_CA.1.F20.C_estpEM_allele_freqs.txt"
# [2] "AS.BF_CA.1.F7.C_estpEM_allele_freqs.txt" 
# [3] "AS.BF_CA.2.F20.C_estpEM_allele_freqs.txt"
# [4] "AS.BF_CA.2.F7.C_estpEM_allele_freqs.txt" 
# [5] "AS.BF_CA.3.F20.C_estpEM_allele_freqs.txt"
# [6] "AS.BF_CA.3.F7.C_estpEM_allele_freqs.txt" 
# [7] "AS.BF_CA.4.F20.C_estpEM_allele_freqs.txt"
# [8] "AS.BF_CA.4.F7.C_estpEM_allele_freqs.txt" 
# [9] "AS.BF_CA.5.F20.C_estpEM_allele_freqs.txt"
#[10] "AS.BF_CA.5.F7.C_estpEM_allele_freqs.txt"
ts<-c(19,6,19,6,19,6,19,6,19,6)
nes<-c(137.7,137.7,156.6,156.6,158.4,158.4,161.4,161.4,150.9,150.9)

p0a<-read.table(aff[fp0[1]],header=FALSE)
p0b<-read.table(aff[fp0[2]],header=FALSE)
p0<-(p0a[,3]+p0b[,3])/2
pvsBFxCA_C<-matrix(NA,nrow=K,ncol=L)
for(k in 1:K){
        pt<-read.table(aff[fpt[k]],header=FALSE)
        pvsBFxCA_C[k,]<-analyNullSimP(p0=p0,p1=pt[,3],t=ts[k],ne=nes[k],lower=TRUE)
}

# get rid of MAF < 0.01 SNPs, too noisy, and non-LG
com<-which(p0 > 0.01 & p0 < 0.99 & is.na(snps[,1])==FALSE)
par(mfrow=c(3,2))
for(k in 1:6){
        plot(-log10(pvsBFxCA_C[k,com]),pch=19)
}

## run PicMin
## generate correlation matrix for null, 10,000 is number of reps
n<-K/2
emp_p_null_dat <- t(replicate(10000, PicMin:::GenerateNullData(1.0, n, a=0.3, b=5, genes=length(com))))

# Calculate the order statistics' p-values for each simulation
emp_p_null_dat_unscaled <- t(apply(emp_p_null_dat ,1, PicMin:::orderStatsPValues))

# Use those p-values to construct the correlation matrix
null_pMax_cor_unscaled <- cor( emp_p_null_dat_unscaled )
null_pMax_cor_unscaled

# convert obs to ordered
opvsBFxCA_C<-pvsBFxCA_C[,com]
for(k in 1:K){
	opvsBFxCA_C[k,]<-PicMin:::EmpiricalPs(opvsBFxCA_C[k,],large_i_small_p=FALSE)
}

# run test and save results
numReps<-20000
picpBFxCA_C_early<-rep(NA,length(com))
piceBFxCA_C_early<-rep(NA,length(com))
picpBFxCA_C_late<-rep(NA,length(com))
piceBFxCA_C_late<-rep(NA,length(com))
early<-seq(2,10,2)
late<-seq(1,9,2)
for(i in 1:length(com)){
	cat(i,"\n")
	o<-PicMin:::PicMin(opvsBFxCA_C[early,i],null_pMax_cor_unscaled, numReps = numReps)
	picpBFxCA_C_early[i]<-o$p
	piceBFxCA_C_early[i]<-o$config_est
	o<-PicMin:::PicMin(opvsBFxCA_C[late,i],null_pMax_cor_unscaled, numReps = numReps)
	picpBFxCA_C_late[i]<-o$p
	piceBFxCA_C_late[i]<-o$config_est
}

## BF X CA on L
fp0<-c(41,73)
fpt<-grep(pattern="AS.BF_CA.[12345].F[0123456789]+.L",x=aff,fixed=FALSE)
K<-length(fpt)
aff[fpt]
# [1] "AS.BF_CA.1.F20.L_estpEM_allele_freqs.txt"
# [2] "AS.BF_CA.1.F7.L_estpEM_allele_freqs.txt" 
# [3] "AS.BF_CA.2.F20.L_estpEM_allele_freqs.txt"
# [4] "AS.BF_CA.2.F7.L_estpEM_allele_freqs.txt" 
# [5] "AS.BF_CA.3.F20.L_estpEM_allele_freqs.txt"
# [6] "AS.BF_CA.3.F7.L_estpEM_allele_freqs.txt" 
# [7] "AS.BF_CA.4.F20.L_estpEM_allele_freqs.txt"
# [8] "AS.BF_CA.4.F7.L_estpEM_allele_freqs.txt" 
# [9] "AS.BF_CA.5.F19.L_estpEM_allele_freqs.txt"
#[10] "AS.BF_CA.5.F7.L_estpEM_allele_freqs.txt" 
ts<-c(19,6,19,6,19,6,19,6,19,6)
nes<-c(116.9,116.9,80,80,78.8,78.8,102.3,102.3,91.1,91.1)

p0a<-read.table(aff[fp0[1]],header=FALSE)
p0b<-read.table(aff[fp0[2]],header=FALSE)
p0<-(p0a[,3]+p0b[,3])/2
pvsBFxCA<-matrix(NA,nrow=K,ncol=L)
for(k in 1:K){
        pt<-read.table(aff[fpt[k]],header=FALSE)
        pvsBFxCA[k,]<-analyNullSimP(p0=p0,p1=pt[,3],t=ts[k],ne=nes[k],lower=TRUE)
}

# get rid of MAF < 0.01 SNPs, too noisy, and non-LG
com<-which(p0 > 0.01 & p0 < 0.99 & is.na(snps[,1])==FALSE)
par(mfrow=c(3,2))
for(k in 1:6){
        plot(-log10(pvsBFxCA[k,com]),pch=19)
}

## run PicMin
## generate correlation matrix for null, 10,000 is number of reps
n<-K/2
emp_p_null_dat <- t(replicate(10000, PicMin:::GenerateNullData(1.0, n, a=0.3, b=5, genes=length(com))))

# Calculate the order statistics' p-values for each simulation
emp_p_null_dat_unscaled <- t(apply(emp_p_null_dat ,1, PicMin:::orderStatsPValues))

# Use those p-values to construct the correlation matrix
null_pMax_cor_unscaled <- cor( emp_p_null_dat_unscaled )
null_pMax_cor_unscaled

# convert obs to ordered
opvsBFxCA<-pvsBFxCA[,com]
for(k in 1:K){
	opvsBFxCA[k,]<-PicMin:::EmpiricalPs(opvsBFxCA[k,],large_i_small_p=FALSE)
}

# run test and save results
numReps<-20000
picpBFxCA_early<-rep(NA,length(com))
piceBFxCA_early<-rep(NA,length(com))
picpBFxCA_late<-rep(NA,length(com))
piceBFxCA_late<-rep(NA,length(com))
early<-seq(2,10,2)
late<-seq(1,9,2)
for(i in 1:length(com)){
	cat(i,"\n")
	o<-PicMin:::PicMin(opvsBFxCA[early,i],null_pMax_cor_unscaled, numReps = numReps)
	picpBFxCA_early[i]<-o$p
	piceBFxCA_early[i]<-o$config_est
	o<-PicMin:::PicMin(opvsBFxCA[late,i],null_pMax_cor_unscaled, numReps = numReps)
	picpBFxCA_late[i]<-o$p
	piceBFxCA_late[i]<-o$config_est
}

## BZ X CA on C
fp0<-c(67,73)
fpt<-grep(pattern="AS.BZ_CA.[12345].F[0123456789]+.C",x=aff,fixed=FALSE)
K<-length(fpt)
aff[fpt]
# [1] "AS.BZ_CA.1.F20.C_estpEM_allele_freqs.txt"
# [2] "AS.BZ_CA.1.F7.C_estpEM_allele_freqs.txt" 
# [3] "AS.BZ_CA.2.F20.C_estpEM_allele_freqs.txt"
# [4] "AS.BZ_CA.2.F7.C_estpEM_allele_freqs.txt" 
# [5] "AS.BZ_CA.3.F20.C_estpEM_allele_freqs.txt"
# [6] "AS.BZ_CA.3.F7.C_estpEM_allele_freqs.txt" 
# [7] "AS.BZ_CA.4.F20.C_estpEM_allele_freqs.txt"
# [8] "AS.BZ_CA.4.F7.C_estpEM_allele_freqs.txt" 
# [9] "AS.BZ_CA.5.F20.C_estpEM_allele_freqs.txt"
#[10] "AS.BZ_CA.5.F7.C_estpEM_allele_freqs.txt" 
ts<-c(19,6,19,6,19,6,19,6,19,6)
nes<-c(191.9,191.9,213.9,213.9,212.8,212.8,142.8,142.8,222.4,222.4)

p0a<-read.table(aff[fp0[1]],header=FALSE)
p0b<-read.table(aff[fp0[2]],header=FALSE)
p0<-(p0a[,3]+p0b[,3])/2
pvsBZxCA_C<-matrix(NA,nrow=K,ncol=L)
for(k in 1:K){
        pt<-read.table(aff[fpt[k]],header=FALSE)
        pvsBZxCA_C[k,]<-analyNullSimP(p0=p0,p1=pt[,3],t=ts[k],ne=nes[k],lower=TRUE)
}

# get rid of MAF < 0.01 SNPs, too noisy, and non-LG
com<-which(p0 > 0.01 & p0 < 0.99 & is.na(snps[,1])==FALSE)
par(mfrow=c(3,2))
for(k in 1:6){
        plot(-log10(pvsBZxCA_C[k,com]),pch=19)
}

## run PicMin
## generate correlation matrix for null, 10,000 is number of reps
n<-K/2
emp_p_null_dat <- t(replicate(10000, PicMin:::GenerateNullData(1.0, n, a=0.3, b=5, genes=length(com))))

# Calculate the order statistics' p-values for each simulation
emp_p_null_dat_unscaled <- t(apply(emp_p_null_dat ,1, PicMin:::orderStatsPValues))

# Use those p-values to construct the correlation matrix
null_pMax_cor_unscaled <- cor( emp_p_null_dat_unscaled )
null_pMax_cor_unscaled

# convert obs to ordered
opvsBZxCA_C<-pvsBZxCA_C[,com]
for(k in 1:K){
	opvsBZxCA_C[k,]<-PicMin:::EmpiricalPs(opvsBZxCA_C[k,],large_i_small_p=FALSE)
}

# run test and save results
numReps<-20000
picpBZxCA_C_early<-rep(NA,length(com))
piceBZxCA_C_early<-rep(NA,length(com))
picpBZxCA_C_late<-rep(NA,length(com))
piceBZxCA_C_late<-rep(NA,length(com))
early<-seq(2,10,2)
late<-seq(1,9,2)
for(i in 1:length(com)){
	cat(i,"\n")
	o<-PicMin:::PicMin(opvsBZxCA_C[early,i],null_pMax_cor_unscaled, numReps = numReps)
	picpBZxCA_C_early[i]<-o$p
	piceBZxCA_C_early[i]<-o$config_est
	o<-PicMin:::PicMin(opvsBZxCA_C[late,i],null_pMax_cor_unscaled, numReps = numReps)
	picpBZxCA_C_late[i]<-o$p
	piceBZxCA_C_late[i]<-o$config_est
}

## BZ X CA on L
fp0<-c(67,73)
fpt<-grep(pattern="AS.BZ_CA.[12345].F[0123456789]+.L",x=aff,fixed=FALSE)
K<-length(fpt)
aff[fpt]
# [1] "AS.BZ_CA.1.F19.L_estpEM_allele_freqs.txt"
# [2] "AS.BZ_CA.1.F6.L_estpEM_allele_freqs.txt" 
# [3] "AS.BZ_CA.2.F20.L_estpEM_allele_freqs.txt"
# [4] "AS.BZ_CA.2.F7.L_estpEM_allele_freqs.txt" 
# [5] "AS.BZ_CA.3.F20.L_estpEM_allele_freqs.txt"
# [6] "AS.BZ_CA.3.F7.L_estpEM_allele_freqs.txt" 
# [7] "AS.BZ_CA.4.F20.L_estpEM_allele_freqs.txt"
# [8] "AS.BZ_CA.4.F7.L_estpEM_allele_freqs.txt" 
# [9] "AS.BZ_CA.5.F20.L_estpEM_allele_freqs.txt"
#[10] "AS.BZ_CA.5.F7.L_estpEM_allele_freqs.txt" 
ts<-c(19,6,19,6,19,6,19,6,18,6)
nes<-c(111.0,111.0,86.7,86.7,90.9,90.9,96.5,96.5,89,89)

p0a<-read.table(aff[fp0[1]],header=FALSE)
p0b<-read.table(aff[fp0[2]],header=FALSE)
p0<-(p0a[,3]+p0b[,3])/2
pvsBZxCA<-matrix(NA,nrow=K,ncol=L)
for(k in 1:K){
        pt<-read.table(aff[fpt[k]],header=FALSE)
        pvsBZxCA[k,]<-analyNullSimP(p0=p0,p1=pt[,3],t=ts[k],ne=nes[k],lower=TRUE)
}

# get rid of MAF < 0.01 SNPs, too noisy, and non-LG
com<-which(p0 > 0.01 & p0 < 0.99 & is.na(snps[,1])==FALSE)
par(mfrow=c(3,2))
for(k in 1:6){
        plot(-log10(pvsBZxCA[k,com]),pch=19)
}

## run PicMin
## generate correlation matrix for null, 10,000 is number of reps
n<-K/2
emp_p_null_dat <- t(replicate(10000, PicMin:::GenerateNullData(1.0, n, a=0.3, b=5, genes=length(com))))

# Calculate the order statistics' p-values for each simulation
emp_p_null_dat_unscaled <- t(apply(emp_p_null_dat ,1, PicMin:::orderStatsPValues))

# Use those p-values to construct the correlation matrix
null_pMax_cor_unscaled <- cor( emp_p_null_dat_unscaled )
null_pMax_cor_unscaled

# convert obs to ordered
opvsBZxCA<-pvsBZxCA[,com]
for(k in 1:K){
	opvsBZxCA[k,]<-PicMin:::EmpiricalPs(opvsBZxCA[k,],large_i_small_p=FALSE)
}

# run test and save results
numReps<-20000
picpBZxCA_early<-rep(NA,length(com))
piceBZxCA_early<-rep(NA,length(com))
picpBZxCA_late<-rep(NA,length(com))
piceBZxCA_late<-rep(NA,length(com))
early<-seq(2,10,2)
late<-seq(1,9,2)
for(i in 1:length(com)){
	cat(i,"\n")
	o<-PicMin:::PicMin(opvsBZxCA[early,i],null_pMax_cor_unscaled, numReps = numReps)
	picpBZxCA_early[i]<-o$p
	piceBZxCA_early[i]<-o$config_est
	o<-PicMin:::PicMin(opvsBZxCA[late,i],null_pMax_cor_unscaled, numReps = numReps)
	picpBZxCA_late[i]<-o$p
	piceBZxCA_late[i]<-o$config_est
}
save(list=ls(),file="null.rdat")


ons<-objects()[63:77] ## need to check each time to make sure range is right
sigPic<-rep(NA,length(ons))
for(i in 1:(length(ons))){
	sigPic[i]<-sum(p.adjust(get(ons[i]),method="fdr") <0.05)
}

## number significant after FDR correction
nsig<-data.frame(ons,sigPic)
#                 ons sigPic
#             picpBF    490
#  picpBFxBZ_C_early   1947
#   picpBFxBZ_C_late   1607
#   picpBFxBZ_early   3177
#     picpBFxBZ_late   2141
#  picpBFxCA_C_early   1993
#  picpBFxCA_C_late   1979
#    picpBFxCA_early   3087
#     picpBFxCA_late   2556
#            picpBZ   1807
# picpBZxCA_C_early    690
#  picpBZxCA_C_late    689
#   picpBZxCA_early    689
#    picpBZxCA_late    835
#            picpCA   1915



sets<-c("BF L","BZ L","CA L","BFxBZ L t1","BFxBZ L t2","BFxBZ C t1","BFxBZ C t2","BFxCA L t1","BFxCA L t2","BFxCA C t1","BFxCA C t2","BZxCA L t1","BZxCA L t2","BZxCA C t1","BZxCA C t2")
sigN<-c(490,1807,1915,3177,2141,1947,1607,3087,2556,1993,1979,689,835,690,689)


#BF = "orange"
#BZ = "red3"
#CA = "purple"
#BF/BZ = "orangered"
#BF/CA = "burlywood4"
#BZ/CA = "mediumvioletred"


## summaries by chromosome and picmin P comparisons across sets of lines
chr_mns<-matrix(NA,nrow=15,ncol=10)
comb_picmin_ps<-matrix(NA,nrow=15,ncol=L)

## bf
fp0<-41
p0<-read.table(aff[fp0],header=FALSE)
com<-which(p0[,3] > 0.01 & p0[,3] < 0.99 & is.na(snps[,1])==FALSE)

chr_mns[1,]<-tapply(INDEX=snps[com,1],X=p.adjust(picpBF,method="fdr") < 0.05,mean)
comb_picmin_ps[1,com]<-picpBF

## bz
fp0<-67
p0<-read.table(aff[fp0],header=FALSE)
com<-which(p0[,3] > 0.01 & p0[,3] < 0.99 & is.na(snps[,1])==FALSE)

chr_mns[2,]<-tapply(INDEX=snps[com,1],X=p.adjust(picpBZ,method="fdr") < 0.05,mean)
comb_picmin_ps[2,com]<-picpBZ

## ca
fp0<-73
p0<-read.table(aff[fp0],header=FALSE)
com<-which(p0[,3] > 0.01 & p0[,3] < 0.99 & is.na(snps[,1])==FALSE)

chr_mns[3,]<-tapply(INDEX=snps[com,1],X=p.adjust(picpCA,method="fdr") < 0.05,mean)
comb_picmin_ps[3,com]<-picpCA

## BF/BZ
p0f<-c(41,67)
p0a<-read.table(aff[p0f[1]],header=FALSE)
p0b<-read.table(aff[p0f[2]],header=FALSE)
p0<-(p0a[,3]+p0b[,3])/2

com<-which(p0 > 0.01 & p0 < 0.99 & is.na(snps[,1])==FALSE)
chr_mns[4,]<-tapply(INDEX=snps[com,1],X=p.adjust(picpBFxBZ_early,method="fdr") < 0.05,mean)
chr_mns[5,]<-tapply(INDEX=snps[com,1],X=p.adjust(picpBFxBZ_late,method="fdr") < 0.05,mean)
chr_mns[6,]<-tapply(INDEX=snps[com,1],X=p.adjust(picpBFxBZ_C_early,method="fdr") < 0.05,mean)
chr_mns[7,]<-tapply(INDEX=snps[com,1],X=p.adjust(picpBFxBZ_C_late,method="fdr") < 0.05,mean)

comb_picmin_ps[4,com]<-picpBFxBZ_early
comb_picmin_ps[5,com]<-picpBFxBZ_late
comb_picmin_ps[6,com]<-picpBFxBZ_C_early
comb_picmin_ps[7,com]<-picpBFxBZ_C_late

## BF/CA
p0f<-c(41,73)
p0a<-read.table(aff[p0f[1]],header=FALSE)
p0b<-read.table(aff[p0f[2]],header=FALSE)
p0<-(p0a[,3]+p0b[,3])/2

com<-which(p0 > 0.01 & p0 < 0.99 & is.na(snps[,1])==FALSE)
chr_mns[8,]<-tapply(INDEX=snps[com,1],X=p.adjust(picpBFxCA_early,method="fdr") < 0.05,mean)
chr_mns[9,]<-tapply(INDEX=snps[com,1],X=p.adjust(picpBFxCA_late,method="fdr") < 0.05,mean)
chr_mns[10,]<-tapply(INDEX=snps[com,1],X=p.adjust(picpBFxCA_C_early,method="fdr") < 0.05,mean)
chr_mns[11,]<-tapply(INDEX=snps[com,1],X=p.adjust(picpBFxCA_C_late,method="fdr") < 0.05,mean)

comb_picmin_ps[8,com]<-picpBFxCA_early
comb_picmin_ps[9,com]<-picpBFxCA_late
comb_picmin_ps[10,com]<-picpBFxCA_C_early
comb_picmin_ps[11,com]<-picpBFxCA_C_late

## BZ/CA
p0f<-c(67,73)
p0a<-read.table(aff[p0f[1]],header=FALSE)
p0b<-read.table(aff[p0f[2]],header=FALSE)
p0<-(p0a[,3]+p0b[,3])/2

com<-which(p0 > 0.01 & p0 < 0.99 & is.na(snps[,1])==FALSE)
chr_mns[12,]<-tapply(INDEX=snps[com,1],X=p.adjust(picpBZxCA_early,method="fdr") < 0.05,mean)
chr_mns[13,]<-tapply(INDEX=snps[com,1],X=p.adjust(picpBZxCA_late,method="fdr") < 0.05,mean)
chr_mns[14,]<-tapply(INDEX=snps[com,1],X=p.adjust(picpBZxCA_C_early,method="fdr") < 0.05,mean)
chr_mns[15,]<-tapply(INDEX=snps[com,1],X=p.adjust(picpBZxCA_C_late,method="fdr") < 0.05,mean)

comb_picmin_ps[12,com]<-picpBZxCA_early
comb_picmin_ps[13,com]<-picpBZxCA_late
comb_picmin_ps[14,com]<-picpBZxCA_C_early
comb_picmin_ps[15,com]<-picpBZxCA_C_late

cors<-matrix(NA,nrow=15,ncol=15)
for(a in 1:15){for(b in 1:15){
	o<-cor.test(-log10(comb_picmin_ps[a,]),-log10(comb_picmin_ps[b,]),na.rm=TRUE)
	cors[a,b]<-o$estimate
}}
rownames(cors)<-sets
colnames(cors)<-sets

## comparing top 5%
summary(apply(chr_mns,1,mean))
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##0.01085 0.01498 0.03368 0.03098 0.04178 0.05223 
sigSnps<-vector("list",15)
for(k in 1:15){
	qq<-quantile(comb_picmin_ps[k,],probs=0.05,na.rm=TRUE)
	sigSnps[[k]]<-which(comb_picmin_ps[k,] <= qq)
}
obs5<-matrix(NA,nrow=15,ncol=15)
null5<-matrix(NA,nrow=15,ncol=15)
for(a in 1:14){for(b in (a+1):15){
	o<-mean(sigSnps[[a]] %in% sigSnps[[b]])
	obs5[a,b]<-o
	nullv<-rep(NA,1000)
	nullv[1]<-o
	for(i in 1:999){
		x<-which(is.na(comb_picmin_ps[b,])==FALSE)
		xx<-sample(x,length(x),replace=FALSE)
		nullb<-rep(NA,L)
		nullb[x]<-comb_picmin_ps[b,xx]
		qq<-quantile(comb_picmin_ps[b,],probs=0.05,na.rm=TRUE)
		nullSig<-which(nullb <= qq)
		nullv[i+1]<-mean(sigSnps[[a]] %in% nullSig)
	}
	null5[a,b]<-mean(nullv)
	null5[b,a]<-mean(nullv >= o)
	obs5[b,a]<-obs5[a,b]/mean(nullv)
}}

rownames(null5)<-sets
colnames(null5)<-sets
rownames(obs5)<-sets
colnames(obs5)<-sets

px<-c(21,19)
cs<-c("gray80","black")
pdf("repeatAmong.pdf",width=8,height=8)
par(mar=c(5,5,.5,.5))
plot(1:15,1:15,type='n',axes=FALSE,xlab="",ylab="")
for(i in 1:14){for(j in (i+1):15){
	points(j,i,pch=19,col=cs[(null5[j,i]<0.05)+1],cex=obs5[j,i]*.5)
}}
axis(1,at=1:15,sets,las=2)
text(1:15,1:15,sets,col=c("orange","red3","purple",rep(c("orangered","burlywood4","mediumvioletred"),each=4)))
legend(2,14,c(.5,1,3,5,9,11),pch=19,pt.cex=(0.5*c(.5,1,3,5,9,11)),bty='n')
dev.off()

pdf("repeatWithin.pdf",width=9.5,height=5)
layout(matrix(1:2,nrow=1,ncol=2),heights=5,widths=c(4,5.5))
par(mar=c(4.5,5,2.5,.5))
dotchart(x=sigN,labels=sets,xlab="Number of loci",pch=19,col=c("orange","red3","purple",rep(c("orangered","burlywood4","mediumvioletred"),each=4)),cex.lab=1.3,cex.axis=1.1)
title(main="(a) Repeated evolution",cex.main=1.3)


ub<-max(as.numeric(chr_mns))
cs<-c("orange","red3","purple",rep(c("orangered","burlywood4","mediumvioletred"),each=4))
lll<-rep(1,15)
lll[c(6,7,10,11,14,15)]<-3
px<-rep(19,15)
px[c(4,6,8,10,12,14)]<-21

plot(1:10,rep(NA,10),type='n',ylim=c(0,1.1*ub),xlab="Chromsome",ylab="Proportion of SNPs",cex.lab=1.3,cex.axis=1.1)
for(k in 1:15){
	points(1:10,chr_mns[k,],type='b',lty=lll[k],col=cs[k],pch=px[k])
}
legend(1,.27,c("BF","BZ","CA","BFxBZ","BFxCA","BZxCA"),fill=unique(cs),bty='n')
legend(3,.27,c("lentil","cowpea"),lty=c(1,3),bty='n')
legend(6,.27,c("t1 (early)","t2 (late)"),pch=c(21,19),bty='n')
title(main="(b) Repeated evolution by chromosome",cex.main=1.3)
dev.off()

############################# null model plots####################3
######non-admixed 

chr<-which(is.na(snps[,1])==FALSE)
bnds<-c(1,which(snps[chr[-1],1] != snps[chr[-72583],1]),72583)

sn<-1:L
nn<-1:72583
mids<-tapply(X=nn,INDEX=snps[chr,1],median)

cm<-1.4;ca<-1.1;cl<-1.4
cx<-c(.3,1)
pdf("nullNonAdmix.pdf",width=8,height=9)
par(mfrow=c(3,1))
par(mar=c(4.5,5,2.5,.5))

## BF on lentil
fp0<-41
p0<-read.table(aff[fp0],header=FALSE)
com<-which(p0[,3] > 0.01 & p0[,3] < 0.99 & is.na(snps[,1])==FALSE)
lp<-as.vector(-log10(pvsBF[,com]))
ub<-max(lp[is.finite(lp)])
K<-dim(pvsBF[,com])[1]
cs<-alpha("orange",rev(seq(.15,.8,length.out=K)))
plot(sn[com],-log10(pvsBF[1,com]),type='n',axes=FALSE,xlab="Chromosome",cex.lab=cl,ylim=c(0,ub),ylab="-log10 P-value")
for(i in seq(1,9,2)){
        polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,ub*1.1,ub*1.1),col=alpha("gray",.5),border=NA)
}
sig<-p.adjust(picpBF,method="fdr") < 0.05
#abline(v=com[a],lwd=.03,col="gray20",lty=3)
for(k in 1:K){
	points(sn[com],-log10(pvsBF[k,com]),pch=19,col=cs[k],cex=cx[sig+1])
}
axis(1,at=mids,1:10,cex.axis=ca)
axis(2,cex.axis=ca)
box()
title(main="(a) Genome-wide evidence of selection for BF on lentil",cex.main=cm)

## BZ on lentil
fp0<-67
p0<-read.table(aff[fp0],header=FALSE)
com<-which(p0[,3] > 0.01 & p0[,3] < 0.99 & is.na(snps[,1])==FALSE)
lp<-as.vector(-log10(pvsBZ[,com]))
ub<-max(lp[is.finite(lp)])
K<-dim(pvsBZ[,com])[1]
cs<-alpha("red3",rev(seq(.15,.8,length.out=K)))
plot(sn[com],-log10(pvsBZ[1,com]),type='n',axes=FALSE,xlab="Chromosome",cex.lab=cl,ylim=c(0,ub),ylab="-log10 P-value")
for(i in seq(1,9,2)){
        polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,ub*1.1,ub*1.1),col=alpha("gray",.5),border=NA)
}
sig<-p.adjust(picpBZ,method="fdr") < 0.05
for(k in 1:K){
	points(sn[com],-log10(pvsBZ[k,com]),pch=19,col=cs[k],cex=cx[sig+1])
}
axis(1,at=mids,1:10,cex.axis=ca)
axis(2,cex.axis=ca)
box()
title(main="(b) Genome-wide evidence of selection for BZ on lentil",cex.main=cm)

## CA on lentil
fp0<-73
p0<-read.table(aff[fp0],header=FALSE)
com<-which(p0[,3] > 0.01 & p0[,3] < 0.99 & is.na(snps[,1])==FALSE)
lp<-as.vector(-log10(pvsCA[,com]))
ub<-max(lp[is.finite(lp)])
K<-dim(pvsCA[,com])[1]
cs<-alpha("purple",rev(seq(.15,.8,length.out=K)))
plot(sn[com],-log10(pvsCA[1,com]),type='n',axes=FALSE,xlab="Chromosome",cex.lab=cl,ylim=c(0,ub),ylab="-log10 P-value")
for(i in seq(1,9,2)){
        polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,ub*1.1,ub*1.1),col=alpha("gray",.5),border=NA)
}
sig<-p.adjust(picpCA,method="fdr") < 0.05
for(k in 1:K){
	points(sn[com],-log10(pvsCA[k,com]),pch=19,col=cs[k],cex=cx[sig+1])
}
axis(1,at=mids,1:10,cex.axis=ca)
axis(2,cex.axis=ca)
box()
title(main="(c) Genome-wide evidence of selection for CA on lentil",cex.main=cm)
dev.off()

############### admix L ###################

pdf("nullAdmix.pdf",width=8,height=9)
par(mfrow=c(3,1))
par(mar=c(4.5,5,2.5,.5))

## BF/BZ
p0f<-c(41,67)
p0a<-read.table(aff[p0f[1]],header=FALSE)
p0b<-read.table(aff[p0f[2]],header=FALSE)
p0<-(p0a[,3]+p0b[,3])/2

late<-seq(1,9,2)

com<-which(p0 > 0.01 & p0 < 0.99 & is.na(snps[,1])==FALSE)
lp<-as.vector(-log10(pvsBFxBZ[late,com]))
ub<-max(lp[is.finite(lp)])
K<-dim(pvsBFxBZ[late,com])[1]
cs<-alpha("orangered",rev(seq(.15,.8,length.out=K)))
plot(sn[com],-log10(pvsBFxBZ[late[1],com]),type='n',axes=FALSE,xlab="Chromosome",cex.lab=cl,ylim=c(0,ub),ylab="-log10 P-value")
for(i in seq(1,9,2)){
        polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,ub*1.1,ub*1.1),col=alpha("gray",.5),border=NA)
}
sig<-p.adjust(picpBFxBZ_late,method="fdr") < 0.05
for(k in 1:K){
	points(sn[com],-log10(pvsBFxBZ[late[k],com]),pch=19,col=cs[k],cex=cx[sig+1])
}
axis(1,at=mids,1:10,cex.axis=ca)
axis(2,cex.axis=ca)
box()
title(main="(a) Evidence of selection in BFxBZ on lentil",cex.main=cm)

## BF/CA on lentil
p0f<-c(41,73)
p0a<-read.table(aff[p0f[1]],header=FALSE)
p0b<-read.table(aff[p0f[2]],header=FALSE)
p0<-(p0a[,3]+p0b[,3])/2

late<-seq(1,9,2)

com<-which(p0 > 0.01 & p0 < 0.99 & is.na(snps[,1])==FALSE)
lp<-as.vector(-log10(pvsBFxCA[late,com]))
ub<-max(lp[is.finite(lp)])
K<-dim(pvsBFxCA[late,com])[1]
cs<-alpha("burlywood4",rev(seq(.15,.8,length.out=K)))
plot(sn[com],-log10(pvsBFxCA[late[1],com]),type='n',axes=FALSE,xlab="Chromosome",cex.lab=cl,ylim=c(0,ub),ylab="-log10 P-value")
for(i in seq(1,9,2)){
        polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,ub*1.1,ub*1.1),col=alpha("gray",.5),border=NA)
}
sig<-p.adjust(picpBFxCA_late,method="fdr") < 0.05
for(k in 1:K){
	points(sn[com],-log10(pvsBFxCA[late[k],com]),pch=19,col=cs[k],cex=cx[sig+1])
}
axis(1,at=mids,1:10,cex.axis=ca)
axis(2,cex.axis=ca)
box()
title(main="(b) Evidence of selection in BFxCA on lentil",cex.main=cm)


## BZ/CA on lentil
p0f<-c(67,73)
p0a<-read.table(aff[p0f[1]],header=FALSE)
p0b<-read.table(aff[p0f[2]],header=FALSE)
p0<-(p0a[,3]+p0b[,3])/2

late<-seq(1,9,2)

com<-which(p0 > 0.01 & p0 < 0.99 & is.na(snps[,1])==FALSE)
lp<-as.vector(-log10(pvsBZxCA[late,com]))
ub<-max(lp[is.finite(lp)])
K<-dim(pvsBZxCA[late,com])[1]
cs<-alpha("mediumvioletred",rev(seq(.15,.8,length.out=K)))
plot(sn[com],-log10(pvsBZxCA[late[1],com]),type='n',axes=FALSE,xlab="Chromosome",cex.lab=cl,ylim=c(0,ub),ylab="-log10 P-value")
for(i in seq(1,9,2)){
        polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,ub*1.1,ub*1.1),col=alpha("gray",.5),border=NA)
}
sig<-p.adjust(picpBZxCA_late,method="fdr") < 0.05
for(k in 1:K){
	points(sn[com],-log10(pvsBZxCA[late[k],com]),pch=19,col=cs[k],cex=cx[sig+1])
}
axis(1,at=mids,1:10,cex.axis=ca)
axis(2,cex.axis=ca)
box()
title(main="(c) Evidence of selection in BZxCA on lentil",cex.main=cm)
dev.off()

############### admix L, early ###################

pdf("nullAdmix_early.pdf",width=8,height=9)
par(mfrow=c(3,1))
par(mar=c(4.5,5,2.5,.5))

## BF/BZ
p0f<-c(41,67)
p0a<-read.table(aff[p0f[1]],header=FALSE)
p0b<-read.table(aff[p0f[2]],header=FALSE)
p0<-(p0a[,3]+p0b[,3])/2

early<-seq(2,10,2)

com<-which(p0 > 0.01 & p0 < 0.99 & is.na(snps[,1])==FALSE)
lp<-as.vector(-log10(pvsBFxBZ[early,com]))
ub<-max(lp[is.finite(lp)])
K<-dim(pvsBFxBZ[early,com])[1]
cs<-alpha("orangered",rev(seq(.15,.8,length.out=K)))
plot(sn[com],-log10(pvsBFxBZ[early[1],com]),type='n',axes=FALSE,xlab="Chromosome",cex.lab=cl,ylim=c(0,ub),ylab="-log10 P-value")
for(i in seq(1,9,2)){
        polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,ub*1.1,ub*1.1),col=alpha("gray",.5),border=NA)
}
sig<-p.adjust(picpBFxBZ_early,method="fdr") < 0.05
for(k in 1:K){
	points(sn[com],-log10(pvsBFxBZ[early[k],com]),pch=19,col=cs[k],cex=cx[sig+1])
}
axis(1,at=mids,1:10,cex.axis=ca)
axis(2,cex.axis=ca)
box()
title(main="(a) Evidence of selection in BFxBZ on lentil",cex.main=cm)

## BF/CA on lentil
p0f<-c(41,73)
p0a<-read.table(aff[p0f[1]],header=FALSE)
p0b<-read.table(aff[p0f[2]],header=FALSE)
p0<-(p0a[,3]+p0b[,3])/2

com<-which(p0 > 0.01 & p0 < 0.99 & is.na(snps[,1])==FALSE)
lp<-as.vector(-log10(pvsBFxCA[early,com]))
ub<-max(lp[is.finite(lp)])
K<-dim(pvsBFxCA[early,com])[1]
cs<-alpha("burlywood4",rev(seq(.15,.8,length.out=K)))
plot(sn[com],-log10(pvsBFxCA[early[1],com]),type='n',axes=FALSE,xlab="Chromosome",cex.lab=cl,ylim=c(0,ub),ylab="-log10 P-value")
for(i in seq(1,9,2)){
        polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,ub*1.1,ub*1.1),col=alpha("gray",.5),border=NA)
}
sig<-p.adjust(picpBFxCA_early,method="fdr") < 0.05
for(k in 1:K){
	points(sn[com],-log10(pvsBFxCA[early[k],com]),pch=19,col=cs[k],cex=cx[sig+1])
}
axis(1,at=mids,1:10,cex.axis=ca)
axis(2,cex.axis=ca)
box()
title(main="(b) Evidence of selection in BFxCA on lentil",cex.main=cm)


## BZ/CA on lentil
p0f<-c(67,73)
p0a<-read.table(aff[p0f[1]],header=FALSE)
p0b<-read.table(aff[p0f[2]],header=FALSE)
p0<-(p0a[,3]+p0b[,3])/2

com<-which(p0 > 0.01 & p0 < 0.99 & is.na(snps[,1])==FALSE)
lp<-as.vector(-log10(pvsBZxCA[early,com]))
ub<-max(lp[is.finite(lp)])
K<-dim(pvsBZxCA[early,com])[1]
cs<-alpha("mediumvioletred",rev(seq(.15,.8,length.out=K)))
plot(sn[com],-log10(pvsBZxCA[early[1],com]),type='n',axes=FALSE,xlab="Chromosome",cex.lab=cl,ylim=c(0,ub),ylab="-log10 P-value")
for(i in seq(1,9,2)){
        polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,ub*1.1,ub*1.1),col=alpha("gray",.5),border=NA)
}
sig<-p.adjust(picpBZxCA_early,method="fdr") < 0.05
for(k in 1:K){
	points(sn[com],-log10(pvsBZxCA[early[k],com]),pch=19,col=cs[k],cex=cx[sig+1])
}
axis(1,at=mids,1:10,cex.axis=ca)
axis(2,cex.axis=ca)
box()
title(main="(c) Evidence of selection in BZxCA on lentil",cex.main=cm)
dev.off()

############### admix C ###################

pdf("nullAdmixC.pdf",width=8,height=9)
par(mfrow=c(3,1))
par(mar=c(4.5,5,2.5,.5))

## BF/BZ
p0f<-c(41,67)
p0a<-read.table(aff[p0f[1]],header=FALSE)
p0b<-read.table(aff[p0f[2]],header=FALSE)
p0<-(p0a[,3]+p0b[,3])/2

late<-seq(1,9,2)

com<-which(p0 > 0.01 & p0 < 0.99 & is.na(snps[,1])==FALSE)
lp<-as.vector(-log10(pvsBFxBZ_C[late,com]))
ub<-max(lp[is.finite(lp)])
K<-dim(pvsBFxBZ_C[late,com])[1]
cs<-alpha("orangered",rev(seq(.15,.8,length.out=K)))
plot(sn[com],-log10(pvsBFxBZ_C[late[1],com]),type='n',axes=FALSE,xlab="Chromosome",cex.lab=cl,ylim=c(0,ub),ylab="-log10 P-value")
for(i in seq(1,9,2)){
        polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,ub*1.1,ub*1.1),col=alpha("gray",.5),border=NA)
}
sig<-p.adjust(picpBFxBZ_C_late,method="fdr") < 0.05
for(k in 1:K){
	points(sn[com],-log10(pvsBFxBZ_C[late[k],com]),pch=19,col=cs[k],cex=cx[sig+1])
}
axis(1,at=mids,1:10,cex.axis=ca)
axis(2,cex.axis=ca)
box()
title(main="(a) Evidence of selection in BFxBZ on cowpea",cex.main=cm)

## BF/CA on cowpea
p0f<-c(41,73)
p0a<-read.table(aff[p0f[1]],header=FALSE)
p0b<-read.table(aff[p0f[2]],header=FALSE)
p0<-(p0a[,3]+p0b[,3])/2

late<-seq(1,9,2)

com<-which(p0 > 0.01 & p0 < 0.99 & is.na(snps[,1])==FALSE)
lp<-as.vector(-log10(pvsBFxCA_C[late,com]))
ub<-max(lp[is.finite(lp)])
K<-dim(pvsBFxCA_C[late,com])[1]
cs<-alpha("burlywood4",rev(seq(.15,.8,length.out=K)))
plot(sn[com],-log10(pvsBFxCA_C[late[1],com]),type='n',axes=FALSE,xlab="Chromosome",cex.lab=cl,ylim=c(0,ub),ylab="-log10 P-value")
for(i in seq(1,9,2)){
        polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,ub*1.1,ub*1.1),col=alpha("gray",.5),border=NA)
}
sig<-p.adjust(picpBFxCA_C_late,method="fdr") < 0.05
for(k in 1:K){
	points(sn[com],-log10(pvsBFxCA_C[late[k],com]),pch=19,col=cs[k],cex=cx[sig+1])
}
axis(1,at=mids,1:10,cex.axis=ca)
axis(2,cex.axis=ca)
box()
title(main="(b) Evidence of selection in BFxCA on cowpea",cex.main=cm)


## BZ/CA on cowpea
p0f<-c(67,73)
p0a<-read.table(aff[p0f[1]],header=FALSE)
p0b<-read.table(aff[p0f[2]],header=FALSE)
p0<-(p0a[,3]+p0b[,3])/2

late<-seq(1,9,2)

com<-which(p0 > 0.01 & p0 < 0.99 & is.na(snps[,1])==FALSE)
lp<-as.vector(-log10(pvsBZxCA_C[late,com]))
ub<-max(lp[is.finite(lp)])
K<-dim(pvsBZxCA_C[late,com])[1]
cs<-alpha("mediumvioletred",rev(seq(.15,.8,length.out=K)))
plot(sn[com],-log10(pvsBZxCA_C[late[1],com]),type='n',axes=FALSE,xlab="Chromosome",cex.lab=cl,ylim=c(0,ub),ylab="-log10 P-value")
for(i in seq(1,9,2)){
        polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,ub*1.1,ub*1.1),col=alpha("gray",.5),border=NA)
}
sig<-p.adjust(picpBZxCA_C_late,method="fdr") < 0.05
for(k in 1:K){
	points(sn[com],-log10(pvsBZxCA_C[late[k],com]),pch=19,col=cs[k],cex=cx[sig+1])
}
axis(1,at=mids,1:10,cex.axis=ca)
axis(2,cex.axis=ca)
box()
title(main="(c) Evidence of selection in BZxCA on cowpea",cex.main=cm)
dev.off()

############### admix C early ###################

pdf("nullAdmixC_early.pdf",width=8,height=9)
par(mfrow=c(3,1))
par(mar=c(4.5,5,2.5,.5))

## BF/BZ
p0f<-c(41,67)
p0a<-read.table(aff[p0f[1]],header=FALSE)
p0b<-read.table(aff[p0f[2]],header=FALSE)
p0<-(p0a[,3]+p0b[,3])/2

early<-seq(2,10,2)

com<-which(p0 > 0.01 & p0 < 0.99 & is.na(snps[,1])==FALSE)
lp<-as.vector(-log10(pvsBFxBZ_C[early,com]))
ub<-max(lp[is.finite(lp)])
K<-dim(pvsBFxBZ_C[early,com])[1]
cs<-alpha("orangered",rev(seq(.15,.8,length.out=K)))
plot(sn[com],-log10(pvsBFxBZ_C[early[1],com]),type='n',axes=FALSE,xlab="Chromosome",cex.lab=cl,ylim=c(0,ub),ylab="-log10 P-value")
for(i in seq(1,9,2)){
        polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,ub*1.1,ub*1.1),col=alpha("gray",.5),border=NA)
}
sig<-p.adjust(picpBFxBZ_C_early,method="fdr") < 0.05
for(k in 1:K){
	points(sn[com],-log10(pvsBFxBZ_C[early[k],com]),pch=19,col=cs[k],cex=cx[sig+1])
}
axis(1,at=mids,1:10,cex.axis=ca)
axis(2,cex.axis=ca)
box()
title(main="(a) Evidence of selection in BFxBZ on cowpea",cex.main=cm)

## BF/CA on cowpea
p0f<-c(41,73)
p0a<-read.table(aff[p0f[1]],header=FALSE)
p0b<-read.table(aff[p0f[2]],header=FALSE)
p0<-(p0a[,3]+p0b[,3])/2

com<-which(p0 > 0.01 & p0 < 0.99 & is.na(snps[,1])==FALSE)
lp<-as.vector(-log10(pvsBFxCA_C[early,com]))
ub<-max(lp[is.finite(lp)])
K<-dim(pvsBFxCA_C[early,com])[1]
cs<-alpha("burlywood4",rev(seq(.15,.8,length.out=K)))
plot(sn[com],-log10(pvsBFxCA_C[early[1],com]),type='n',axes=FALSE,xlab="Chromosome",cex.lab=cl,ylim=c(0,ub),ylab="-log10 P-value")
for(i in seq(1,9,2)){
        polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,ub*1.1,ub*1.1),col=alpha("gray",.5),border=NA)
}
sig<-p.adjust(picpBFxCA_C_early,method="fdr") < 0.05
for(k in 1:K){
	points(sn[com],-log10(pvsBFxCA_C[early[k],com]),pch=19,col=cs[k],cex=cx[sig+1])
}
axis(1,at=mids,1:10,cex.axis=ca)
axis(2,cex.axis=ca)
box()
title(main="(b) Evidence of selection in BFxCA on cowpea",cex.main=cm)


## BZ/CA on cowpea
p0f<-c(67,73)
p0a<-read.table(aff[p0f[1]],header=FALSE)
p0b<-read.table(aff[p0f[2]],header=FALSE)
p0<-(p0a[,3]+p0b[,3])/2

com<-which(p0 > 0.01 & p0 < 0.99 & is.na(snps[,1])==FALSE)
lp<-as.vector(-log10(pvsBZxCA_C[early,com]))
ub<-max(lp[is.finite(lp)])
K<-dim(pvsBZxCA_C[early,com])[1]
cs<-alpha("mediumvioletred",rev(seq(.15,.8,length.out=K)))
plot(sn[com],-log10(pvsBZxCA_C[early[1],com]),type='n',axes=FALSE,xlab="Chromosome",cex.lab=cl,ylim=c(0,ub),ylab="-log10 P-value")
for(i in seq(1,9,2)){
        polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,ub*1.1,ub*1.1),col=alpha("gray",.5),border=NA)
}
sig<-p.adjust(picpBZxCA_C_early,method="fdr") < 0.05
for(k in 1:K){
	points(sn[com],-log10(pvsBZxCA_C[early[k],com]),pch=19,col=cs[k],cex=cx[sig+1])
}
axis(1,at=mids,1:10,cex.axis=ca)
axis(2,cex.axis=ca)
box()
title(main="(c) Evidence of selection in BZxCA on cowpea",cex.main=cm)
dev.off()
