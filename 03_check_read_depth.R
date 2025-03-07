#This is a small script to take a look at the distribution of read depths in genomic data. 
#The purpose of checking read depth is to set an appropriate maximum read depth for the second filtering step to remove probable orthologous sites (i.e. sites that have double or triple (etc.) the average read depth because two separate orthologous regions of the genome mapped to a single reference). 

bob <- read.table("cmac_admix_read_depth.txt", stringsAsFactors = FALSE, header = FALSE)
bob <- as.numeric(bob[,1])

hist(bob, breaks = 50)

summary(bob)

#C. mac. admxiture stats: 

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 3072    5214    8147   12394   15361  333328 

#Hayden's ringlet data for comparison: 

#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#771    1488    1986    2582    2967   20750 

SD <- sd(bob)
#1811.531 for ringlets
#11979.01 for cmac admix

mean(bob) + 3*(SD)
#48330.97

#Thus, will be using around 48,000 as max read depth for the C. maculatus admixture data. 
