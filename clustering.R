library(kernlab)
library(ggplot2)
require(plot3D)
source("coSequenceCPP.R")
source("gapkernel.R")

set.seed(6)
N= 50
reuters <- read.table("reuters.txt.gz", header=T)
reuters <- reuters[reuters$Topic == "crude" | reuters$Topic == "grain" | reuters$Topic == "coffee",]
reuters <- reuters[sample(1:nrow(reuters),N),]
reuters$Content <- as.character(reuters$Content)  
reuters$Topic <- factor(reuters$Topic)            
levels(reuters$Topic)
length(reuters$Topic)
table(reuters$Topic)

getClusters <- function( kernelMatrix, k, bool ){
  sc <- kkmeans(K, reuters$Topic, centers = k, nstart = 10)
  if(bool==TRUE){
    vec <-vector(length = N)
    for(x in 1:N){ vec[x] = sc[x]}
    list(table(vec, reuters$Topic), withinss(sc), sum(withinss(sc)) )
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
writeInFile<- function(text,R){
  write(text,file = 'ResultsClustering.txt', append=TRUE)
  write.table( R[[1]],file = 'ResultsClustering.txt', append=TRUE ,sep = ",")
  write(R[[2]],file = 'ResultsClustering.txt', append=TRUE)
  write(R[[3]],file = 'ResultsClustering.txt', append=TRUE)
}

write( paste("Experiment with ",as.character(N)," texts"),file = 'ResultsClustering.txt', append=TRUE)
write.table(table(reuters$Topic),file = 'ResultsClustering.txt', append=TRUE ,sep = ",")
clusters = 3
set.seed(20)
n = bestSubSeqLength("spectrum", reuters$Content, clusters)
k <- stringdot("spectrum",length=n)
K <- kernelMatrix(k, reuters$Content)
R<-getClusters(K,clusters,T)
writeInFile("SPECTRUM KERNEL",R)

k <- stringdot("CONSTANT")
K <- kernelMatrix(k, reuters$Content)
R<-getClusters(K,clusters,T)
writeInFile("CONSTANT KERNEL",R)

n = bestSubSeqLength("boundrange", reuters$Content, clusters)
k <- stringdot("boundrange",length=n)
K <- kernelMatrix(k, reuters$Content)
R<-getClusters(K,clusters,T)
writeInFile("BOUNDRANGE KERNEL",R)

lamdas = seq(1.2,2,0.2)
M <- bestDecayFactor("exponential",reuters$Content,clusters, lamdas)
k <- stringdot("exponential",length= M[2], lambda = lamdas[M[1]] )
K <- kernelMatrix(k, reuters$Content)
R<-getClusters(K,clusters,T)
writeInFile("EXPONENTIAL KERNEL",R)


bestPossibleClustering <- function(lambdas, n , data ,k ){
  M <- matrix(nrow = length(lambdas), ncol = n-1)
  for(i in 1:length(lambdas)) {
    M[i,] =sapply(seq(2,n), function(x){ kp = makeCppKernel(lambdas[i],x) 
                                         K <- kernelMatrix(kp, data)
                                         getClusters(K,k,FALSE) })
  }
  which(M == min(M), arr.ind = TRUE) 
}

lambdas = seq(0.1,0.9,0.1)
fc = bestPossibleClustering(lambdas, 4, reuters$Content,clusters)
k <- makeCppKernel(lambdas[fc[1]], fc[2])
K <- kernelMatrix(k, reuters$Content)
writeInFile("GAP KERNEL",R)

k <- new("kernel", .Data=coSequenceKernelCPP, kpar=list())
K <- kernelMatrix(k ,reuters$Content)
R <-getClusters(K,clusters,T)
writeInFile("COSEQUENCE KERNEL",R)

