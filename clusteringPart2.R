## This part asks for the discovering of the
## correct number of clusters.
library(kernlab)
library(ggplot2)
require(plot3D)
soucre("clustering.R")
set.seed(6)

findBestK <- function(K, limit,name, kReal){
  ss <-vector(length= limit)
  for( i in 1: limit) ss[i] <- getClusters(K,i)
  k = which.min(ss)
  jpeg('name')
  plot(seq(1,limit), ss, 
       main = paste(name,"-kernel ", "k optimal = ",as.character(k),"k real = ",as.character(kReal) ),
       xlab = "different k", ylab="sum of ss")
  dev.off()
}

N = 30
reuters <- read.table("reuters.txt.gz", header=T)
kReal = 10
topics <- sample(unique(reuters[1,]),k)

reuters <- reuters[topics,]
reuters <- reuters[sample(1:nrow(reuters),N),]

reuters$Content <- as.character(reuters$Content)  
reuters$Topic <- factor(reuters$Topic)            

table(reuters$Topic)

k <- new("kernel", .Data=coSequenceKernelCPP, kpar=list())
K <- kernelMatrix(k ,reuters$Content)
findBestK(K,50,"CONNECT",kReal ) 

