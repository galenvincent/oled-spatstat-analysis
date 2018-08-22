# Test file for the HPC
print('Lets get things going...')

library(parallel)

print('parallel library sucessfully uploaded')

#library(spatstat)
#print('spatstat library sucessfully uploaded')

print(detectCores())

a <- c(1:20)

write.table(a, file="output.csv",sep=",",row.names=F)
