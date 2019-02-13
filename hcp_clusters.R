# Creation of HCP clusters to check the K function model 
library(spatstat)
library(rapt)
library(data.table)

###################################################################

# HCP cluster with Poisson underneath
# define volume domain 
hcp.vol <- box3(c(0,60),c(0,60),c(0,60))
cluster <- hcpcluster(6,3,0.5,0.03,hcp.vol,'poisson')
plot3d.pp3(cluster[[1]])
cluster[[2]]

b <- G3est(cluster[[1]])
plot(b)

a <- K3est(cluster[[1]],rmax = 11,nrval = 300, correction = "translation")
plot(a,xlim=c(0,11))
a.df <- as.data.frame(cbind(a$r,a$trans))

# HCP cluster with RCP underneath
hcp.vol <- box3(c(0,60),c(0,60),c(0,60))
cluster <- hcpcluster(6,3,0.5,0.03,hcp.vol,'rcp',c("~/Research/point_patterns/Final/FinalConfig1","~/Research/point_patterns/Final/system1"))
#plot3d.pp3(cluster[[1]])
cluster[[2]]

a <- K3est(cluster[[1]],rmax = 11,nrval = 300, correction = "translation")
#plot(a,xlim=c(0,11))
a.df <- as.data.frame(cbind(a$r,a$trans))

fwrite(a.df,"hcpclust_rcp_csep_r6-R3-s105-s2003.csv")

# get the 50,000 RCP RRL ready to plot in mathematica 
toSub <- fread('~/Research/K_cluster_series/cubetoSub.csv',drop=1)
env.r <- fread('~/Research/K_cluster_series/cube.csv',select=2)
sim.r <- c(0, 11/19, 22/19, 33/19, 44/19, 55/19, 66/19, 77/19, 88/19, 99/19, 
110/19, 121/19, 132/19, 143/19, 154/19, 165/19, 176/19, 187/19, 
198/19, 11)

toSub.new <- rep(0,20)

for(i in 1:20){
  toSub.new[i] <- toSub$x[which(abs(env.r - sim.r[i]) == min(abs(env.r - sim.r[i])))]
}
fwrite(as.data.frame(toSub.new),"analyticModel_toSub.csv")
