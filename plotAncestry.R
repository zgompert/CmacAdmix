## plot popanc results
library(scales)


pf<-list.files(path="output_summary_files/",full.names=TRUE,pattern="d_1000000_n_15_popQ")
K<-5

cl<-1.4;cm<-1.4;ca<-1.1

pdf("ancestryFreqt2.pdf",width=8,height=9)
par(mfrow=c(3,1))
par(mar=c(4.5,5,2.5,.5))


## BFxBZ t2 (late)
lf<-c(1,6,10,14,18)
cf<-c(2,5,9,13,17)

snps<-scan("input_popanc_files/scaf_BFxBZ.txt")
snps[snps==1]<-NA
chr<-as.numeric(as.factor(snps))
a<-which(is.na(chr)==FALSE)
bnds<-c(1,which(chr[a][-1] != chr[a][-length(a)]),length(a))
mids<-tapply(X=1:length(a),INDEX=chr[a],median)

cs<-"orangered"


plot(1:length(a),rep(.5,length(a)),ylim=c(0,1),axes=FALSE,type='n',xlab="Chromosome",ylab="BF ancestry",cex.lab=cl)

for(i in seq(1,9,2)){
        polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,1.1,1.1),col=alpha("gray",.5),border=NA)
}
lancP<-matrix(NA,nrow=5,ncol=length(a))
cancP<-matrix(NA,nrow=5,ncol=length(a))
mnsBFxBZ_L<-matrix(NA,nrow=5,ncol=10)
mnsBFxBZ_C<-matrix(NA,nrow=5,ncol=10)
for(k in 1:5){
	lanc<-read.table(pf[lf[k]],header=TRUE,sep=",")
	canc<-read.table(pf[cf[k]],header=TRUE,sep=",")
	lancP[k,]<-lanc$median[a]
	cancP[k,]<-canc$median[a]
	mnsBFxBZ_L[k,]<-tapply(X=lancP[k,],INDEX=chr[a],median)
	mnsBFxBZ_C[k,]<-tapply(X=cancP[k,],INDEX=chr[a],median)
}
abline(h=.5,lty=2)
lines(apply(cancP,2,median),lty=1,col=alpha(cs,.5))
lines(apply(lancP,2,median),lty=1,col=cs)

legend(1,1.05,c(paste("lentil; mean = ",round(mean(lancP),2),sep=""),paste("cowpea; mean = ",round(mean(cancP),2),sep="")),lwd=2,col=c(cs,alpha(cs,.5)),bty='n')

axis(1,at=mids,1:10,cex.axis=ca)
axis(2,at=c(0,.5,1),cex.axis=ca)
title(main="(a) Genome-wide ancestry frequencies for BFxBZ",cex.main=cm)
box()


## BFxCA t2 (late)
lf<-c(2,6,10,14,17)+20
cf<-c(1,5,9,13,18)+20

snps<-scan("input_popanc_files/scaf_BFxCA.txt")
snps[snps==1]<-NA
chr<-as.numeric(as.factor(snps))
a<-which(is.na(chr)==FALSE)
bnds<-c(1,which(chr[a][-1] != chr[a][-length(a)]),length(a))
mids<-tapply(X=1:length(a),INDEX=chr[a],median)

cs<-"burlywood4"


plot(1:length(a),rep(.5,length(a)),ylim=c(0,1),axes=FALSE,type='n',xlab="Chromosome",ylab="BF ancestry",cex.lab=cl)

for(i in seq(1,9,2)){
        polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,1.1,1.1),col=alpha("gray",.5),border=NA)
}
lancP<-matrix(NA,nrow=5,ncol=length(a))
cancP<-matrix(NA,nrow=5,ncol=length(a))
mnsBFxCA_L<-matrix(NA,nrow=5,ncol=10)
mnsBFxCA_C<-matrix(NA,nrow=5,ncol=10)
for(k in 1:5){
	lanc<-read.table(pf[lf[k]],header=TRUE,sep=",")
	canc<-read.table(pf[cf[k]],header=TRUE,sep=",")
	lancP[k,]<-lanc$median[a]
	cancP[k,]<-canc$median[a]
	mnsBFxCA_L[k,]<-tapply(X=lancP[k,],INDEX=chr[a],median)
	mnsBFxCA_C[k,]<-tapply(X=cancP[k,],INDEX=chr[a],median)
}
abline(h=.5,lty=2)
lines(apply(cancP,2,median),lty=1,col=alpha(cs,.5))
lines(apply(lancP,2,median),lty=1,col=cs)

legend(1,1.05,c(paste("lentil; mean = ",round(mean(lancP),2),sep=""),paste("cowpea; mean = ",round(mean(cancP),2),sep="")),lwd=2,col=c(cs,alpha(cs,.5)),bty='n')

axis(1,at=mids,1:10,cex.axis=ca)
axis(2,at=c(0,.5,1),cex.axis=ca)
title(main="(b) Genome-wide ancestry frequencies for BFxCA",cex.main=cm)
box()


## BZxCA t2 (late)
lf<-c(1,6,10,14,18)+40
cf<-c(2,5,9,13,17)+40

snps<-scan("input_popanc_files/scaf_BZxCA.txt")
snps[snps==1]<-NA
chr<-as.numeric(as.factor(snps))
a<-which(is.na(chr)==FALSE)
bnds<-c(1,which(chr[a][-1] != chr[a][-length(a)]),length(a))
mids<-tapply(X=1:length(a),INDEX=chr[a],median)

cs<-"mediumvioletred"


plot(1:length(a),rep(.5,length(a)),ylim=c(0,1),axes=FALSE,type='n',xlab="Chromosome",ylab="BZ ancestry",cex.lab=cl)

for(i in seq(1,9,2)){
        polygon(c(bnds[i],bnds[i+1],bnds[i+1],bnds[i]),c(-.1,-.1,1.1,1.1),col=alpha("gray",.5),border=NA)
}
lancP<-matrix(NA,nrow=5,ncol=length(a))
cancP<-matrix(NA,nrow=5,ncol=length(a))
mnsBZxCA_L<-matrix(NA,nrow=5,ncol=10)
mnsBZxCA_C<-matrix(NA,nrow=5,ncol=10)
for(k in 1:5){
	lanc<-read.table(pf[lf[k]],header=TRUE,sep=",")
	canc<-read.table(pf[cf[k]],header=TRUE,sep=",")
	lancP[k,]<-lanc$median[a]
	cancP[k,]<-canc$median[a]
	mnsBZxCA_L[k,]<-tapply(X=lancP[k,],INDEX=chr[a],median)
	mnsBZxCA_C[k,]<-tapply(X=cancP[k,],INDEX=chr[a],median)
}
abline(h=.5,lty=2)
lines(apply(cancP,2,median),lty=1,col=alpha(cs,.5))
lines(apply(lancP,2,median),lty=1,col=cs)

legend(1,1.05,c(paste("lentil; mean = ",round(mean(lancP),2),sep=""),paste("cowpea; mean = ",round(mean(cancP),2),sep="")),lwd=2,col=c(cs,alpha(cs,.5)),bty='n')

axis(1,at=mids,1:10,cex.axis=ca)
axis(2,at=c(0,.5,1),cex.axis=ca)
title(main="(c) Genome-wide ancestry frequencies for BZxCA",cex.main=cm)
box()
dev.off()


pdf("chromAncestryFreqt2.pdf",width=6.5,height=7.3)
par(mfrow=c(3,1))
par(mar=c(4.5,5,2.5,.5))

cs<-"orangered"
xx<-barplot(rbind(mnsBFxBZ_C,mnsBFxBZ_L),beside=TRUE,col=rep(c(alpha(cs,.5),cs),each=5),ylab="BF ancestry",xlab="Chromosome",ylim=c(0,.8),cex.lab=cl)
axis(1,at=(xx[5,]+xx[6,])/2,1:10,cex.axis=ca)
abline(h=.5,lty=2)
box()
legend(0,.82,c("cowpea","lentil"),fill=c(alpha(cs,.5),cs),bty='n')

title(main="(a) Ancestry by chromosome for BFxBZ",cex.main=cm)




cs<-"burlywood4"
xx<-barplot(rbind(mnsBFxCA_C,mnsBFxCA_L),beside=TRUE,col=rep(c(alpha(cs,.5),cs),each=5),ylab="BF ancestry",xlab="Chromosome",ylim=c(0,.8),cex.lab=cl)
axis(1,at=(xx[5,]+xx[6,])/2,1:10,cex.axis=ca)
abline(h=.5,lty=2)
box()
legend(0,.82,c("cowpea","lentil"),fill=c(alpha(cs,.5),cs),bty='n')

title(main="(b) Ancestry by chromosome for BFxCA",cex.main=cm)

cs<-"mediumvioletred"
xx<-barplot(rbind(mnsBZxCA_C,mnsBZxCA_L),beside=TRUE,col=rep(c(alpha(cs,.5),cs),each=5),ylab="BZ ancestry",xlab="Chromosome",ylim=c(0,.8),cex.lab=cl)
axis(1,at=(xx[5,]+xx[6,])/2,1:10,cex.axis=ca)
abline(h=.5,lty=2)
box()
legend(0,.82,c("cowpea","lentil"),fill=c(alpha(cs,.5),cs),bty='n')

title(main="(c) Ancestry by chromosome for BZxCA",cex.main=cm)
dev.off()
