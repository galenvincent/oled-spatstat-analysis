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

under.nums <- seq(2,102,1)
under.nums[length(under.nums)] <- 1
over.nums <- seq(1,101,1)

hcp.vol <- box3(c(0,60),c(0,60),c(0,60))

envout <- matrix(NA, nrow = 300, ncol = 102)

for(i in 1:101){
  print(i)
  fp <- c(paste('~/Research/point_patterns/Final/FinalConfig',toString(under.nums[i]),sep=""),paste('~/Research/point_patterns/Final/system',toString(under.nums[i]),sep=""))
  cluster <- hcpcluster(6,3,0.5,0.03,hcp.vol,'rcp',fp)
  #plot3d.pp3(cluster[[1]])
  #cluster[[2]]
  
  a <- K3est(cluster[[1]],rmax = 11,nrval = 300, correction = "translation")
  if(i == 1){
    envout[,i] <- a$r
  }
  envout[,i+1] <- a$trans
}

envout.vals <- envout[,2:102]
envout.sorted <- t(apply(envout.vals, 1, sort))

perc <- 0.95
envibig <- ceiling((1-((1-perc)/2))*100)
envismall <- floor(((1-perc)/2)*100)

env.toplot <- cbind(envout[,1],envout.sorted[,envismall],envout.sorted[,envibig])

View(env.toplot)

a.df <- as.data.frame(env.toplot)

fwrite(a.df,"hcpclust_rcp_env_csep_r6-R3-s105-s2003.csv")




# get the 50,000 RCP RRL ready to plot in mathematica 
toSub <- fread('~/Research/numeric_k_model/cubetoSub_r11.csv',drop=1)
env.r <- fread('~/Research/numeric_k_model/cube_r11.csv',select=2)
sim.r <- seq(0,11,len = 30)

toSub.new <- rep(0,30)
toSub.r <- rep(0,30)

for(i in 1:30){
  toSub.r[i] <- env.r$V1[which(abs(env.r - sim.r[i]) == min(abs(env.r - sim.r[i])))]
  toSub.new[i] <- toSub$x[which(abs(env.r - sim.r[i]) == min(abs(env.r - sim.r[i])))]
}
fwrite(as.data.frame(toSub.new),"~/Research/numeric_k_model/analyticModel_toSub.csv")
fwrite(as.data.frame(toSub.r),"~/Research/numeric_k_model/analyticModel_toSub_r.csv")
