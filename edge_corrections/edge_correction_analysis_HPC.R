# Comparing border vs trans & iso edge corrections
library(rapt)
library(parallel)

# Upload an RCP pattern
r1 <-0.0252429 #1.2 sigma
rcp_1_upload <- read.table("C:/Users/galen/OneDrive/Documents/Research/point_patterns/Final/FinalConfig1",sep = " ",col.names = c("x","y","z","type"))
rcp_1 <- scaleRCP(createSpat(rcp_1_upload[,c("x","y","z")]),newRadius = 0.5,oldRadius = r1)

rcp <- stitch(rcp_1)
# Run random relabeling envelopes with border correction on the full sample up to about 5nm radius
maxR <- 10
bordRRL <- pK3est(0.06,rcp,1000,rmax=maxR,correction="bord",anom=TRUE)
write.table(bordRRL[[1]], file = "C:/Users/galen/Desktop/bord.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")

# Run random relabeling envelopes with trans on smaller inner box sample
inner.size <- domain(rcp)$xrange[2]-maxR*2
inner.box <- subSquare(rcp,rep(inner.size,3))

transRRL <- pK3est(0.06,inner.box,1000,rmax=maxR,correction="trans",anom=TRUE,toSub=bordRRL[[2]])
write.table(transRRL[[1]], file = "C:/Users/galen/Desktop/trans.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")

isoRRL <- pK3est(0.06,inner.box,1000,rmax=maxR,correction="iso",anom=TRUE,toSub=bordRRL[[2]])
write.table(isoRRL[[1]], file = "C:/Users/galen/Desktop/iso.csv",row.names=FALSE, na="",col.names=FALSE, sep=",")