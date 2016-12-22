library(kernlab)
library(tm)
source("connect.R")
source("gapkernel.R")


set.seed(6)
N= 10 # number of data point
reuters <- read.table("reuters.txt.gz", header=TRUE)
reuters <- reuters[reuters$Topic == "crude" | reuters$Topic == "grain" | reuters$Topic == "coffee",]
reuters <- reuters[sample(1:nrow(reuters),N),]
reuters$Content <- as.character(reuters$Content)  
reuters$Topic <- factor(reuters$Topic)

docs <- Corpus(VectorSource(reuters$Content))
docs <- tm_map(docs, removePunctuation)
docs <- tm_map(docs, stripWhitespace)
docs <- tm_map(docs,removeNumbers)
reuters$Content <- sapply(docs, function(x){x$content})

levels(reuters$Topic)
length(reuters$Topic)
table(reuters$Topic)

getClusters <- function( kernelMatrix, k, bool ){
  sc <- kkmeans(kernelMatrix, reuters$Topic, centers = k, nstart = 10)
  if(bool==TRUE){
    vec <-vector(length = N)
    for(x in 1:N){ vec[x] = sc[x]}
    list(table(vec, reuters$Topic), withinss(sc), sum(withinss(sc)) )
  }
  if(bool == FALSE){
    sum (withinss(sc))
  }
}
impactK <- function( M ,s ,l ,k , data, ... ){
  vec <- sapply(seq(2,s), function(x){ 
                    kp <- stringdot(M, length = x, lambda = l)
                    K <- kernelMatrix( kp, data )
                    getClusters(K, k, FALSE) })

}

minLength<- function( TextSet ){
  min( sapply( TextSet, function(x){ nchar(x) } ) )
}

bestSubSeqLength <- function( Kernel, data, clusters, ... ){
  values <- impactK(Kernel, minLength(data), 1 , clusters, data)
  which.min(values)
}

bestDecayFactor <- function( Kernel, Data, cluster, range ){
  M <- matrix(nrow = length(range), ncol = minLength(Data)-1 )
  for(i in 1:length(range)) {
    M[i,] <- impactK(Kernel, minLength(Data), range[i], clusters, Data)}
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

##Calculate the clusteres using the spectrum kernel
n = bestSubSeqLength("spectrum", reuters$Content, clusters)
k <- stringdot("spectrum",length=n)
K <- kernelMatrix(k, reuters$Content)
R <- getClusters(K,clusters,TRUE)
writeInFile("SPECTRUM KERNEL",R)
##Calculate the clusteres using the constant kernel
k <- stringdot("constant",length=2)
K <- kernelMatrix(k, reuters$Content)
R<-getClusters(K,clusters,TRUE)
writeInFile("CONSTANT KERNEL",R)
##Calculate the clusters using the bound range kernel
n = bestSubSeqLength("boundrange", reuters$Content, clusters)
k <- stringdot("boundrange",length=n)
K <- kernelMatrix(k, reuters$Content)
R<-getClusters(K,clusters,T)
writeInFile("BOUNDRANGE KERNEL",R)
##Calculate the cluster using the exponential kernel
lamdas = seq(1.2,2,0.2)
M <- bestDecayFactor("exponential",reuters$Content,clusters, lamdas)
k <- stringdot("exponential",length= M[2], lambda = lamdas[M[1]] )
K <- kernelMatrix(k, reuters$Content)
R<-getClusters(K,clusters,T)
writeInFile("EXPONENTIAL KERNEL",R)

##Function to find good parameters for the gap kernel
bestPossibleClustering <- function(lambdas, n , data ,k ){
  M <- matrix(nrow = length(lambdas), ncol = n-1)
  for(i in 1:length(lambdas)) {
    M[i,] =sapply(seq(2,n), function(x){ kp = makeCppKernel(lambdas[i],x) 
                                         K <- kernelMatrix(kp, data)
                                         getClusters(K,k,FALSE) })
  }
  which(M == min(M), arr.ind = TRUE) 
}
##Launch the gap kernel
lambdas = seq(0.1,0.9,0.1)
fc = bestPossibleClustering(lambdas, 6, reuters$Content,clusters)
k <- makeCppKernel(lambdas[fc[1]], fc[2])
K <- kernelMatrix(k, reuters$Content)
writeInFile("GAP KERNEL",R)
##Clustering using the connect kernel
k <- new("kernel", .Data=connect, kpar=list())
K <- kernelMatrix(k ,reuters$Content)
R <-getClusters(K,clusters,TRUE)
writeInFile("CONNECT",R)

