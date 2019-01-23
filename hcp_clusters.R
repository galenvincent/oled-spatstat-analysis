# Creation of HCP clusters to check the K function model 
library(spatstat)
library(rapt)
library(data.table)

###################################################################

# define volume domain 
hcp.vol <- box3(c(0,60),c(0,60),c(0,60))
cluster <- hcpcluster(6,3,0.5,0.03,hcp.vol)
plot(cluster[[1]])
cluster[[2]]

a <- K3est(cluster[[1]],rmax = 11,nrval = 300, correction = "translation")
plot(a,xlim=c(0,11))
a.df <- as.data.frame(cbind(a$r,a$trans))
fwrite(a.df,"hcpclust_csep_r6-R3-s105-s2003.csv")

