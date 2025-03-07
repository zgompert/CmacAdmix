#This is a script for A. Springer's Callosobruchus maculatus admixture adaptation project.
#This takes raw count data for the number of beetles produced each 20 days 
#and plots population growth and cumulative jar productivity

#####################################################################################
## Load libraries
#####################################################################################
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())



#######################################################################################
## READ IN DATA 
#######################################################################################

data <- read.csv("Lentil_pop_growth_data.csv", stringsAsFactors = FALSE)

rep.colors <- c("firebrick3", "darkorange", "gold", "chartreuse", "green4", 
                "blue", "deepskyblue", "darkviolet", "magenta3", "hotpink", "peru")




#######################################################################################
## CALCULATE CUMULATIVE COUNTS
#######################################################################################


#calculate the cumulative productivity of each jar:
cumulative_count <- numeric(length = nrow(data))
for(i in 1:nrow(data)){
  cumulative_count[i] <- sum(data[which(data$Population == data$Population[i] & 
                                          data$Replicate == data$Replicate[i] &
                                          data$Total_number_days_post_host_shift <= 
                                          data$Total_number_days_post_host_shift[i]),16])
  #print(data[which(data$Population == data$Population[i] & 
  #data$Replicate == data$Replicate[i] &
  #data$Total_number_days_post_host_shift <= 
  #data$Total_number_days_post_host_shift[i]),16])
}


#######################################################################################
## PLOT CUMULATIVE COUNT DATA
#######################################################################################



pdf("cmac_admix_cum_pop_counts.pdf", width = 9, height = 6)

##MAKE TWO PLOTS IN ONE IMAGE:
par(mfrow = c(2, 3))

##set figure margins (bottom, left, top, right):
par(mar = c(5, 6, 5, 2))

#BURKINA FASO CUMULATIVE:
plot(NULL, 
     xlim = c(0,400), ylim = c(0,20000),
     xlab = "No. Days Post Host Shift", ylab = "Total Beetles Emerged",
     main = "BF Cumulative Pop. vs. Time")
for(i in 6:11){
  points(x = data[which(data$Replicate == i & data$Population == "BF"),10], 
         y = cumulative_count[which(data$Replicate == i & data$Population == "BF")],
         col = "orange", pch = 15-6+i)
}

#BRAZIL CUMULATIVE:
plot(NULL, 
     xlim = c(0,400), ylim = c(0,20000),
     xlab = "No. Days Post Host Shift", ylab = "Total Beetles Emerged",
     main = "BZ Cumulative Pop. vs. Time")
for(i in 6:11){
  points(x = data[which(data$Replicate == i & data$Population == "BZ"),10], 
         y = cumulative_count[which(data$Replicate == i & data$Population == "BZ")],
         col = "red3", pch = 15-6+i)
}

#CALIFORNIA CUMULATIVE:
plot(NULL, 
     xlim = c(0,400), ylim = c(0,20000),
     xlab = "No. Days Post Host Shift", ylab = "Total Beetles Emerged",
     main = "CA Cumulative Pop. vs. Time")
for(i in 6:11){
  points(x = data[which(data$Replicate == i & data$Population == "CA"),10], 
         y = cumulative_count[which(data$Replicate == i & data$Population == "CA")],
         col = "purple", pch = 15-6+i)
}

#BURKINA FASO/BRAZIL CUMULATIVE:
plot(NULL, 
     xlim = c(0,400), ylim = c(0,20000),
     xlab = "No. Days Post Host Shift", ylab = "Total Beetles Emerged",
     main = "BF/BZ Cumulative Pop. vs. Time")
for(i in 6:11){
  points(x = data[which(data$Replicate == i & data$Population == "BF/BZ"),10], 
         y = cumulative_count[which(data$Replicate == i & data$Population == "BF/BZ")],
         col = "orangered", pch = 15-6+i)
}

#BURKINA FASO/CALIFORNIA CUMULATIVE:
plot(NULL, 
     xlim = c(0,400), ylim = c(0,20000),
     xlab = "No. Days Post Host Shift", ylab = "Total Beetles Emerged",
     main = "BF/CA Cumulative Pop. vs. Time")
for(i in 6:11){
  points(x = data[which(data$Replicate == i & data$Population == "BF/CA"),10], 
         y = cumulative_count[which(data$Replicate == i & data$Population == "BF/CA")],
         col = "burlywood4", pch = 15-6+i)
}

#BRAZIL/CALIFORNIA CUMULATIVE:
plot(NULL, 
     xlim = c(0,400), ylim = c(0,20000),
     xlab = "No. Days Post Host Shift", ylab = "Total Beetles Emerged",
     main = "BZ/CA Cumulative Pop. vs. Time")
for(i in 6:11){
  points(x = data[which(data$Replicate == i & data$Population == "BZ/CA"),10], 
         y = cumulative_count[which(data$Replicate == i & data$Population == "BZ/CA")],
         col = "mediumvioletred", pch = 15-6+i)
}


dev.off()


#######################################################################################
## CALCULATE CUMULATIVE COUNTS
#######################################################################################



##  TRY A HEIRARCHICAL POLYNOMIAL REGRESSION FOR MY CUMULATIVE DATA: 


##need to get data ready for model

## N is a scalar representing the total number of data points I have across reps, pops, and dates 
## count is the cumulative number of beetles collected from each jar at each date
## days is a vector of length N representing the date each population count was colleccted on DATA
## combo is a vector of length N representing the pop x rep combo of each data point
## pop is a vector of length N representing the pop of each data point
## Npops is the number of populations (in this case, 6
## Ncombos is a scalar representing the total number of pop x rep combos (in this case, 36)
## poptime is a vector of length N representing the pop x time combination of each data point
## days_poptime is a vector of length Npoptime representing the days post host-shift for each unique pop x time combo

## pull out just the population data for reps 6:11 
pop <- as.numeric(as.factor(data$Population[which(data$Replicate >5)]))
## create a combo vector with pop x rep combo information for each data point
combo <- as.numeric(as.factor(paste(data$Population[which(data$Replicate > 5)],
                                    data$Replicate[which(data$Replicate > 5)])))
##create poptime vector with pop x time combo information for each data point
poptime <- as.numeric(as.factor(paste(data$Population[which(data$Replicate >5)], 
                                      data$Total_number_days_post_host_shift[which(data$Replicate > 5)])))
## only unique pop x time combo information
unique_poptime <- unique(paste(data$Population[which(data$Replicate >5)], 
                               data$Total_number_days_post_host_shift[which(data$Replicate > 5)]))
## the days post host shift information for each pop x time combo
## this was done to ensure that I didn't lose data 
## regarding the exact date each count was collected 
## since some were 60 days, others 67, etc. 
days_poptime <- character(length = length(unique(poptime)))
for(i in 1:length(unique(poptime))){
  days_poptime[i] <- strsplit(unique_poptime, split = " ")[[i]][2]
}
days_poptime <- as.numeric(days_poptime)


########################################################################################
## Zach's code for the analysis starts here
########################################################################################

## give SD 1, but don't center
stnd<-function(x=NA){
	x<-x/sd(x,na.rm=TRUE)
	return(x)
}

## lines 6+
lns<-which(data$Replicate > 5)
## matrix of indicator variables for population effectss
wghts<-matrix(0,nrow=length(lns),ncol=6)
wghts[grep(x=data$Population[lns],pattern="BF"),1]<-1
wghts[grep(x=data$Population[lns],pattern="BZ"),2]<-1
wghts[grep(x=data$Population[lns],pattern="CA"),3]<-1
for(i in 1:length(lns)){
	wghts[i,1:3]<-wghts[i,1:3]/sum(wghts[i,1:3]) ## normalize ind pop effects
}
wghts[grep(x=data$Population[lns],pattern="BF/BZ"),4]<-1
wghts[grep(x=data$Population[lns],pattern="BF/CA"),5]<-1
wghts[grep(x=data$Population[lns],pattern="BZ/CA"),6]<-1

csum<-data$Count[lns]
pops<-unique(data$Population)
reps<-c(6:11)
for(i in 1:6){
	for(j in 1:6){
		a<-which(data$Population[lns]==pops[i] & data$Replicate[lns]==reps[j])
		csum[a]<-cumsum(csum[a])
	}}

reps<-as.numeric(as.factor(data$Replicate[lns]-5))## replicate number 

model_data<-list(Y=csum,day=stnd(data$Total_number_days_post_host_shift[lns]),day2=stnd(data$Total_number_days_post_host_shift[lns])^2,pops=wghts,N=length(lns),Np=6,Nr=6,reps=reps)

n_chains<-5

initf<-function(Np=3,chain_id=1){
        list(Bpop=runif(Np,-10,10),BpopD=runif(Np,-2,2),BpopD2=runif(Np,-.2,.2),alpha=chain_id)
}
init_ll <- lapply(1:n_chains, function(id) initf(chain_id = id))

fit<-stan("popsize.stan",data=model_data,chains=n_chains,iter=3000,warmup=1500)
## diagnostics look great

## compute r2
sig<-extract(fit,"sigma")
vart<-var(model_data$Y) ## total variance
r2<-1-(sig[[1]]^2/vart)
quantile(r2,probs=c(.5,.025,.975))
#      50%      2.5%     97.5%
#0.9259324 0.9163227 0.9342282
## quite high, note without line effects it was in the upper .7s

## focus on core effects, betas, can look at alphas of course, but don't really "care" about them, reps are different
## order is BF, BZ, CA, BF/BZ, BF/CA, BZ/CA
bd<-extract(fit,"BpopD")
bd2<-extract(fit,"BpopD2")

pops<-c("BF","BZ","CA","BF/BZ","BF/CA","BZ/CA")
bdest<-apply(bd[[1]],2,quantile,probs=c(.5,.025,.975))
bd2est<-apply(bd2[[1]],2,quantile,probs=c(.5,.025,.975))
colnames(bdest)<-pops
colnames(bd2est)<-pops
round(bdest,1)
#            BF     BZ    CA  BF/BZ  BF/CA  BZ/CA
#  50%   -145.2  -76.6 519.7  -61.5  -21.0   54.4
#  2.5%  -291.5 -229.6 373.1 -230.4 -183.9 -102.5
#  97.5%    3.0   68.5 666.5  102.0  144.4  215.0
round(bd2est,1)
#           BF    BZ    CA BF/BZ BF/CA BZ/CA
#  50%    81.3 245.7 394.5 501.8 757.7 210.8
#  2.5%   30.9 195.1 345.4 431.1 692.6 153.2
#  97.5% 131.1 298.6 443.5 573.0 821.2 267.9

## bd is essentially telling you about the rate of population growth at 
## day 0. Point estimates range from negative to positive, but CIs are only positive
## for CA. Thus, for most, it isn't clear whether they would persist initially.

## bd2 is positive if the rate of growth increases, negative if decreases. Positive would mean adaptation is happening.
## All are positive. Importantly, the hybrid effects add on to the population effects, and thus all hybrids
## are credibly increasing more than non-hybrids

save(list=ls(),file="popsize.rdat")
