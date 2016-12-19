## This part asks for the discovering of the
## correct number of clusters.
library(kernlab)
library(ggplot2)
require(plot3D)

set.seed(6)

findBestK <- function(K, limit,name, kReal){
  ss <-vector(length = limit)
  for( i in seq(2,limit+1)) {ss[i-1] <- sum(withinss(kkmeans(K, centers = i, nstart = 10)))}
  k = which.min(ss)
  jpeg('name')
  plot(seq(1,limit), ss, 
       main = paste(name,"-kernel ", "k optimal = ",as.character(k),"k real = ",as.character(kReal) ),
       xlab = "different k", ylab="sum of ss")
  dev.off()
}

N = 30

reuters <- read.table("reuters.txt.gz", header=T)
reuters$Content <- as.character(reuters$Content)  
reuters$Topic <- factor(reuters$Topic)

cat = 15
topics <- factor( sample(unique(reuters$Topic),cat) )
reuters <- reuters[ sapply(reuters[,1], function (x) x %in% topics) , ]
reuters <- reuters[sample(1:nrow(reuters),N),]
kReal = length(unique(reuters[,1]))

set.seed(20)

limit = 10
k <- stringdot("exponential",length= 6, lambda = 2 )
K <- kernelMatrix(k, reuters$Content)
findBestK(K, limit,"EXPONENTIAL",kReal ) 

k <- stringdot("CONSTANT")
K <- kernelMatrix(k, reuters$Content) 
findBestK(K,limit,"CONSTANT",kReal ) 

k <- new("kernel", .Data=coSequenceKernelCPP, kpar=list())
K <- kernelMatrix(k ,reuters$Content)
findBestK(K,limit,"CONNECT",kReal ) 


k <- makeCppKernel(lambdas[fc[1]], fc[2])
K <- kernelMatrix(k, reuters$Content)
findBestK(K,limit,"GAP",kReal ) 


