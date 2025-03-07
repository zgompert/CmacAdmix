## make raw allele frequency change plots
library(scales)
library(RColorBrewer)
wavg<-function(x,sz=50){
	tl<-length(x)
	wav<-rep(NA,tl-sz)
	#wav<-matrix(NA,nrow=tl-sz,ncol=2)
	for(i in 1:(tl-sz)){
		wav[i,]<-mean(x[i:(i+sz-1)])
		#wav[i,]<-quantile(x[i:(i+sz-1)],probs=c(.1,.9))
	}
	wav
}
wl<-50

#cs<-brewer.pal(5,"Pastel2")

## files
aff<-list.files(path="output_allele_freq_files",pattern="AS.*",full.names=TRUE)

#BF = "orange"
#BZ = "red3"
#CA = "purple"
#BF/BZ = "orangered"
#BF/CA = "burlywood4"
#BZ/CA = "mediumvioletred"

## metadata
L<-79079
N<-length(aff)

snps<-read.table("../15_drift_sim/ChSNPs.txt")

chr<-which(is.na(snps[,1])==FALSE)
#    1     2     3     4     5     6     7     8     9    10 
# 8198  3609  5242  8727  9287  8093  7178  2884  7506 11859 
length(chr)
#[1] 72583
bnds<-c(1,which(snps[chr[-1],1] != snps[chr[-72583],1]),72583)


nn<-1:72583
mids<-tapply(X=nn,INDEX=snps[chr,1],median)

cm<-1.4;ca<-1.1;cl<-1.4

pdf("changeNonAdmix.pdf",width=8,height=9)
par(mfrow=c(3,1))
par(mar=c(4.5,5,2.5,.5))

## BF
p0f<-41
pLf<-42:45 ## no line 5
K<-length(pLf)

p0<-read.table(aff[p0f],header=FALSE)
pL<-matrix(NA,nrow=K,ncol=L)
dp<-pL
for(k in 1:K){
	pL[k,]<-read.table(aff[pLf[k]],header=FALSE)[,3]
	dp[k,]<-abs(p0[,3]-pL[k,])
}

cs<-alpha("orange",rev(seq(.15,.8,length.out=K)))

mn<-mean(as.numeric(dp[,chr]))
#[1] 0.09431991

plot(dp[1,chr],type='n',ylim=c(0,1),axes=FALSE,xlab="Chromosome",ylab="Allele frequency change",cex.lab=cl)
for(i in seq(1,9,2)){
	polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,1.1,1.1),col=alpha("gray",.5),border=NA)
}
for(k in 1:K){
#	wav<-wavg(dp[k,chr],wl)
#	a<-which(dp[k,chr] > quantile(dp[k,chr],probs=.999))
	a<-which(dp[k,chr] > 0.03)
	points(nn[a],dp[k,chr][a],pch=19,col=cs[k])
#	lb<-min(nn)+wl/2;ub<-max(nn)-wl/2
#	polygon(c(seq(lb,ub,length.out=dim(wav)[1]),rev(seq(lb,ub,length.out=dim(wav)[1]))),c(wav[,1],wav[,2]),border=NA,col=alpha(cs[k],.3))
	#lines(nn[-c(1:wl)]+wl/2,wav,col=cs[k],lwd=2)
}
text(mids[10],1,paste("mean = ",round(mn,3),sep=""),cex=ca)
axis(1,at=mids,1:10,cex.axis=ca)
axis(2,at=seq(0,1,.25),cex.axis=ca)
box()
title(main="(a) Genome-wide change on lentil for BF",cex.main=cm)

## BZ
p0f<-67
pLf<-68:72
K<-length(pLf)

p0<-read.table(aff[p0f],header=FALSE)
pL<-matrix(NA,nrow=K,ncol=L)
dp<-pL
for(k in 1:K){
	pL[k,]<-read.table(aff[pLf[k]],header=FALSE)[,3]
	dp[k,]<-abs(p0[,3]-pL[k,])
}

cs<-alpha("red3",rev(seq(.15,.8,length.out=K)))

mn<-mean(as.numeric(dp[,chr]))
#[1] 0.06392623

plot(dp[1,chr],type='n',ylim=c(0,1),axes=FALSE,xlab="Chromosome",ylab="Allele frequency change",cex.lab=cl)
for(i in seq(1,9,2)){
	polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,1.1,1.1),col=alpha("gray",.5),border=NA)
}
for(k in 1:K){
	#wav<-wavg(dp[k,chr],wl)
	#a<-which(dp[k,chr] > quantile(dp[k,chr],probs=.999))
	a<-which(dp[k,chr] > 0.03)
	points(nn[a],dp[k,chr][a],pch=19,col=cs[k])
	#lb<-min(nn)+wl/2;ub<-max(nn)-wl/2
	#polygon(c(seq(lb,ub,length.out=dim(wav)[1]),rev(seq(lb,ub,length.out=dim(wav)[1]))),c(wav[,1],wav[,2]),border=NA,col=alpha(cs[k],.3))
	#lines(nn[-c(1:wl)]+wl/2,wav,col=cs[k],lwd=2)
}
text(mids[10],1,paste("mean = ",round(mn,3),sep=""),cex=ca)
axis(1,at=mids,1:10,cex.axis=ca)
axis(2,at=seq(0,1,.25),cex.axis=ca)
box()
title(main="(b) Genome-wide change on lentil for BZ",cex.main=cm)

## CA
p0f<-73
pLf<-74:78
K<-length(pLf)

p0<-read.table(aff[p0f],header=FALSE)
pL<-matrix(NA,nrow=K,ncol=L)
dp<-pL
for(k in 1:K){
	pL[k,]<-read.table(aff[pLf[k]],header=FALSE)[,3]
	dp[k,]<-abs(p0[,3]-pL[k,])
}

cs<-alpha("purple",rev(seq(.15,.8,length.out=K)))

mn<-mean(as.numeric(dp[,chr]))

plot(dp[1,chr],type='n',ylim=c(0,1),axes=FALSE,xlab="Chromosome",ylab="Allele frequency change",cex.lab=cl)
for(i in seq(1,9,2)){
	polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,1.1,1.1),col=alpha("gray",.5),border=NA)
}
for(k in 1:K){
	#wav<-wavg(dp[k,chr],wl)
	a<-which(dp[k,chr] > 0.03)
	#a<-which(dp[k,chr] > quantile(dp[k,chr],probs=.999))
	points(nn[a],dp[k,chr][a],pch=19,col=cs[k])
	#lb<-min(nn)+wl/2;ub<-max(nn)-wl/2
	#polygon(c(seq(lb,ub,length.out=dim(wav)[1]),rev(seq(lb,ub,length.out=dim(wav)[1]))),c(wav[,1],wav[,2]),border=NA,col=alpha(cs[k],.3))
	#lines(nn[-c(1:wl)]+wl/2,wav,col=cs[k],lwd=2)
}
text(mids[10],1,paste("mean = ",round(mn,3),sep=""),cex=ca)
axis(1,at=mids,1:10,cex.axis=ca)
axis(2,at=seq(0,1,.25),cex.axis=ca)
box()
title(main="(c) Genome-wide change on lentil for CA",cex.main=cm)
dev.off()

################### admix L ############

pdf("changeAdmixLearly.pdf",width=8,height=9)
par(mfrow=c(3,1))
par(mar=c(4.5,5,2.5,.5))

## BF/BZ
p0f<-c(41,67)
pLf<-c(3,8,12,16,20)
K<-length(pLf)

p0a<-read.table(aff[p0f[1]],header=FALSE)
p0b<-read.table(aff[p0f[2]],header=FALSE)
p0<-(p0a[,3]+p0b[,3])/2
pL<-matrix(NA,nrow=K,ncol=L)
dp<-pL
for(k in 1:K){
	pL[k,]<-read.table(aff[pLf[k]],header=FALSE)[,3]
	dp[k,]<-abs(p0-pL[k,])
}

cs<-alpha("orangered",rev(seq(.15,.8,length.out=K)))

mn<-mean(as.numeric(dp[,chr]))

plot(dp[1,chr],type='n',ylim=c(0,1),axes=FALSE,xlab="Chromosome",ylab="Allele frequency change",cex.lab=cl)
for(i in seq(1,9,2)){
	polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,1.1,1.1),col=alpha("gray",.5),border=NA)
}
for(k in 1:K){
	a<-which(dp[k,chr] > 0.03)
	points(nn[a],dp[k,chr][a],pch=19,col=cs[k])
}
text(mids[10],1,paste("mean = ",round(mn,3),sep=""),cex=ca)
axis(1,at=mids,1:10,cex.axis=ca)
axis(2,at=seq(0,1,.25),cex.axis=ca)
box()
title(main="(a) Genome-wide change on lentil for BFxBZ",cex.main=cm)

## BF/CA
p0f<-c(41,73)
pLf<-c(24,28,32,36,40)
K<-length(pLf)

p0a<-read.table(aff[p0f[1]],header=FALSE)
p0b<-read.table(aff[p0f[2]],header=FALSE)
p0<-(p0a[,3]+p0b[,3])/2

pL<-matrix(NA,nrow=K,ncol=L)
dp<-pL
for(k in 1:K){
	pL[k,]<-read.table(aff[pLf[k]],header=FALSE)[,3]
	dp[k,]<-abs(p0-pL[k,])
}

cs<-alpha("burlywood4",rev(seq(.15,.8,length.out=K)))

mn<-mean(as.numeric(dp[,chr]))

plot(dp[1,chr],type='n',ylim=c(0,1),axes=FALSE,xlab="Chromosome",ylab="Allele frequency change",cex.lab=cl)
for(i in seq(1,9,2)){
	polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,1.1,1.1),col=alpha("gray",.5),border=NA)
}
for(k in 1:K){
	a<-which(dp[k,chr] > 0.03)
	points(nn[a],dp[k,chr][a],pch=19,col=cs[k])
}
text(mids[10],1,paste("mean = ",round(mn,3),sep=""),cex=ca)
axis(1,at=mids,1:10,cex.axis=ca)
axis(2,at=seq(0,1,.25),cex.axis=ca)
box()
title(main="(b) Genome-wide change on lentil for BFxCA",cex.main=cm)

## BZ/CA
p0f<-c(67,73)
pLf<-c(49,54,58,62,66)
K<-length(pLf)

p0a<-read.table(aff[p0f[1]],header=FALSE)
p0b<-read.table(aff[p0f[2]],header=FALSE)
p0<-(p0a[,3]+p0b[,3])/2

pL<-matrix(NA,nrow=K,ncol=L)
dp<-pL
for(k in 1:K){
	pL[k,]<-read.table(aff[pLf[k]],header=FALSE)[,3]
	dp[k,]<-abs(p0-pL[k,])
}

cs<-alpha("mediumvioletred",rev(seq(.15,.8,length.out=K)))

mn<-mean(as.numeric(dp[,chr]))

plot(dp[1,chr],type='n',ylim=c(0,1),axes=FALSE,xlab="Chromosome",ylab="Allele frequency change",cex.lab=cl)
for(i in seq(1,9,2)){
	polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,1.1,1.1),col=alpha("gray",.5),border=NA)
}
for(k in 1:K){
	a<-which(dp[k,chr] > 0.03)
	points(nn[a],dp[k,chr][a],pch=19,col=cs[k])
}
text(mids[10],1,paste("mean = ",round(mn,3),sep=""),cex=ca)
axis(1,at=mids,1:10,cex.axis=ca)
axis(2,at=seq(0,1,.25),cex.axis=ca)
box()
title(main="(c) Genome-wide change on lentil for BZxCA",cex.main=cm)
dev.off()

################### admix C ############

pdf("changeAdmixCearly.pdf",width=8,height=9)
par(mfrow=c(3,1))
par(mar=c(4.5,5,2.5,.5))

## BF/BZ
p0f<-c(41,67)
pLf<-c(4,7,11,15,19)
K<-length(pLf)

p0a<-read.table(aff[p0f[1]],header=FALSE)
p0b<-read.table(aff[p0f[2]],header=FALSE)
p0<-(p0a[,3]+p0b[,3])/2
pL<-matrix(NA,nrow=K,ncol=L)
dp<-pL
for(k in 1:K){
	pL[k,]<-read.table(aff[pLf[k]],header=FALSE)[,3]
	dp[k,]<-abs(p0-pL[k,])
}

cs<-alpha("orangered",rev(seq(.15,.8,length.out=K)))

mn<-mean(as.numeric(dp[,chr]))

plot(dp[1,chr],type='n',ylim=c(0,1),axes=FALSE,xlab="Chromosome",ylab="Allele frequency change",cex.lab=cl)
for(i in seq(1,9,2)){
	polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,1.1,1.1),col=alpha("gray",.5),border=NA)
}
for(k in 1:K){
	a<-which(dp[k,chr] > 0.03)
	points(nn[a],dp[k,chr][a],pch=19,col=cs[k])
}
text(mids[10],1,paste("mean = ",round(mn,3),sep=""),cex=ca)
axis(1,at=mids,1:10,cex.axis=ca)
axis(2,at=seq(0,1,.25),cex.axis=ca)
box()
title(main="(a) Genome-wide change on cowpea for BFxBZ",cex.main=cm)

## BF/CA
p0f<-c(41,73)
pLf<-c(23,27,31,35,39)
K<-length(pLf)

p0a<-read.table(aff[p0f[1]],header=FALSE)
p0b<-read.table(aff[p0f[2]],header=FALSE)
p0<-(p0a[,3]+p0b[,3])/2

pL<-matrix(NA,nrow=K,ncol=L)
dp<-pL
for(k in 1:K){
	pL[k,]<-read.table(aff[pLf[k]],header=FALSE)[,3]
	dp[k,]<-abs(p0-pL[k,])
}

cs<-alpha("burlywood4",rev(seq(.15,.8,length.out=K)))

mn<-mean(as.numeric(dp[,chr]))

plot(dp[1,chr],type='n',ylim=c(0,1),axes=FALSE,xlab="Chromosome",ylab="Allele frequency change",cex.lab=cl)
for(i in seq(1,9,2)){
	polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,1.1,1.1),col=alpha("gray",.5),border=NA)
}
for(k in 1:K){
	a<-which(dp[k,chr] > 0.03)
	points(nn[a],dp[k,chr][a],pch=19,col=cs[k])
}
text(mids[10],1,paste("mean = ",round(mn,3),sep=""),cex=ca)
axis(1,at=mids,1:10,cex.axis=ca)
axis(2,at=seq(0,1,.25),cex.axis=ca)
box()
title(main="(b) Genome-wide change on cowpea for BFxCA",cex.main=cm)

## BZ/CA
p0f<-c(67,73)
pLf<-c(50,53,57,61,65)
K<-length(pLf)

p0a<-read.table(aff[p0f[1]],header=FALSE)
p0b<-read.table(aff[p0f[2]],header=FALSE)
p0<-(p0a[,3]+p0b[,3])/2

pL<-matrix(NA,nrow=K,ncol=L)
dp<-pL
for(k in 1:K){
	pL[k,]<-read.table(aff[pLf[k]],header=FALSE)[,3]
	dp[k,]<-abs(p0-pL[k,])
}

cs<-alpha("mediumvioletred",rev(seq(.15,.8,length.out=K)))

mn<-mean(as.numeric(dp[,chr]))

plot(dp[1,chr],type='n',ylim=c(0,1),axes=FALSE,xlab="Chromosome",ylab="Allele frequency change",cex.lab=cl)
for(i in seq(1,9,2)){
	polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,1.1,1.1),col=alpha("gray",.5),border=NA)
}
for(k in 1:K){
	a<-which(dp[k,chr] > 0.03)
	points(nn[a],dp[k,chr][a],pch=19,col=cs[k])
}
text(mids[10],1,paste("mean = ",round(mn,3),sep=""),cex=ca)
axis(1,at=mids,1:10,cex.axis=ca)
axis(2,at=seq(0,1,.25),cex.axis=ca)
box()
title(main="(c) Genome-wide change on cowpea for BZxCA",cex.main=cm)
dev.off()


