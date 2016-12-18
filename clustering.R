library(kernlab)
library(ggplot2)
require(plot3D)
source("coSequenceCPP.R")
source("gapkernel.R")

set.seed(6)
n = 30
reuters <- read.table("reuters.txt.gz", header=T)
reuters <- reuters[reuters$Topic == "crude" | reuters$Topic == "grain" | reuters$Topic == "coffee",]
reuters <- reuters[sample(1:nrow(reuters),n),]
reuters$Content <- as.character(reuters$Content)  
reuters$Topic <- factor(reuters$Topic)            
levels(reuters$Topic)
length(reuters$Topic)
table(reuters$Topic)

getClusters <- function( kernelMatrix, k, bool ){
  sc <- kkmeans(K, reuters$Topic, centers = k, nstart = 10)
  if(bool==TRUE){
    for(x in 1:30){ vec[x] = sc[x]}
    table(vec, reuters$Topic) 
    #print ( withinss(sc) )
  }else {
    sum (withinss(sc))
  }
}
impactK <- function( M,s,l,k,data, ... ){
  vec <- sapply(seq(2,s), function(x){ 
                    kp <- stringdot(M, length = x, lambda = l)
                    K <- kernelMatrix(kp,data)
                    getClusters(K,k,FALSE) })

}
minLength<- function( TextSet ){
  min( sapply( TextSet, function(x){ nchar(x) } ) )
}
bestSubSeqLength <- function( Kernel,Data,clusters, ... ){
  values <- impactK(Kernel, minLength(Data),1, clusters, Data)
  plot(seq(1,length(values)),values, 
       main = paste(Kernel,"-kernel"),xlab = "length n", ylab="sum of ss")
  
  which.min(values)
}
bestDecayFactor <- function( Kernel, Data, cluster, range ){
  M <- matrix(nrow = length(range),ncol = minLength(Data) -1 )
  for(i in 1:length(range)) {
    M[i,] <- impactK(Kernel, minLength(Data), range[i], clusters, Data)}
  persp3D(z = M, theta = 120)
  which(M == min(M), arr.ind = TRUE) 
}

write( paste("Experiment with ",as.character(n)," texts"),file = 'ResultsClustering.txt', append=TRUE)
clusters = 3
set.seed(20)
n = bestSubSeqLength("spectrum", reuters$Content, clusters)
k <- stringdot("spectrum",length=5)
K <- kernelMatrix(k, reuters$Content)
R<-getClusters(K,clusters,T)
write("SPECTRUM KERNEL",file = 'ResultsClustering.txt', append=TRUE)
write.table( R,file = 'ResultsClustering.txt', append=TRUE ,sep = ",")

n = bestSubSeqLength("boundrange", reuters$Content, clusters)
k <- stringdot("boundrange",length=n)
K <- kernelMatrix(k, reuters$Content)
R<-getClusters(K,clusters,T)
write("BOUNDRANGE KERNEL",file = 'ResultsClustering.txt', append=TRUE)
write.table( R,file = 'ResultsClustering.txt', append=TRUE ,sep = ",")

lamdas = seq(1.2,2,0.2)
M <- bestDecayFactor("exponential",reuters$Content,clusters, lamdas)
k <- stringdot("exponential",length= M[2], lambda = lamdas[M[1]] )
K <- kernelMatrix(k, reuters$Content)
R<-getClusters(K,clusters,T)
write("EXPONENTIAL KERNEL",file = 'ResultsClustering.txt', append=TRUE)
write.table( R,file = 'ResultsClustering.txt', append=TRUE ,sep = ",")


bestPossibleClustering <- function(lambdas, n , data ,k ){
  M <- matrix(nrow = length(lambdas), ncol = n-1)
  for(i in 1:length(range)) {
    M[i,] =sapply(seq(2,n), function(x){ kp = makeCppKernel(lambdas[i],x) 
                                         K <- kernelMatrix(kp, data)
                                         getClusters(K,k,FALSE) })
  }
  which(M == min(M), arr.ind = TRUE) 
}

lambdas = seq(0.1,0.9,0.1)
fc = bestPossibleClustering(lambdas, 10, reuters$Content,clusters)
k <- makeCppKernel(lambdas[fc[1]], fc[2])
K <- kernelMatrix(k3, reuters$Content)
R<-getClusters(K,clusters,T)
write("GAP KERNEL",file = 'ResultsClustering.txt', append=TRUE)
write.table( R,file = 'ResultsClustering.txt', append=TRUE ,sep = ",")

k <- new("kernel", .Data=coSequenceKernelCPP, kpar=list())
K <- kernelMatrix(k ,reuters$Content)
R <-getClusters(K,clusters,T)
write("COSEQUENCE KERNEL",file = 'ResultsClustering.txt', append=TRUE)
write.table( R,file = 'ResultsClustering.txt', append=TRUE ,sep = ",")

