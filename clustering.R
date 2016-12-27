library(kernlab)
library(tm)
source("connect.R")
source("gapkernel.R")

#Reading and sampling of the data
N= 50 # number of data point
reuters <- read.table("reuters.txt.gz", header=TRUE)
reuters <- reuters[reuters$Topic == "crude" | reuters$Topic == "grain" | reuters$Topic == "coffee",]
reuters <- reuters[sample(1:nrow(reuters),N),]
reuters$Content <- as.character(reuters$Content)  
reuters$Topic <- factor(reuters$Topic)
## Preprocessing of the data
docs <- Corpus(VectorSource(reuters$Content))
docs <- tm_map(docs, removePunctuation)
docs <- tm_map(docs, stripWhitespace)
docs <- tm_map(docs,removeNumbers)
reuters$Content <- sapply(docs, function(x){x$content})

## function that cluster the data with matrix and params
cluster <- function(kernelMatrix, k, output){
  set.seed(20)
  sc <- kkmeans(kernelMatrix, k, nstart = 10)
  if(output){
    vec <-vector(length = N)
    for(x in 1:N){ vec[x] = sc[x] }
    return( list( table( vec, reuters$Topic ), withinss(sc), sum( withinss(sc) ) ) )
  }
  if(!output){
    return (sum (withinss(sc)) )
  }
}
## shortest text length in the sample
maxLength <- function( TextSet ){
  min( sapply( TextSet, function(x){ nchar(x) } ) )
}
## try a possible string length
bestLength <- function(kernelName, k, dataSet){
  max <- 15#maxLength(dataSet)
  scores <- sapply( seq(2,max), function(x) { param <- stringdot(kernelName, length = x)
                                              kernelMatrix <- kernelMatrix(param,dataSet)
                                              possibleError <- tryCatch(
                                                sc <- cluster(kernelMatrix, k, FALSE),
                                                error=function(e) return(NA)
                                              )
                                              if(!inherits(possibleError, "error")){
                                                return (sc)
                                              }
                                              
                                             })
  return( which.min(scores) + 1 )
}

writeInFile<- function(text,R){
  suppressWarnings(write(text,file = 'ResultsClustering.txt', append=TRUE))
  suppressWarnings(write.table( R[[1]],file = 'ResultsClustering.txt', append=TRUE ,sep = ","))
  suppressWarnings(write(R[[2]],file = 'ResultsClustering.txt', append=TRUE))
  suppressWarnings(write(R[[3]],file = 'ResultsClustering.txt', append=TRUE))
}

suppressWarnings(write( paste("Experiment with ",as.character(N)," texts"),file = 'ResultsClustering.txt', append=TRUE))
suppressWarnings(write.table(table(reuters$Topic),file = 'ResultsClustering.txt', append=TRUE ,sep = ","))

clusters = 3
##Calculate the clusters using spectrum kernel
n <- bestLength("spectrum",clusters,reuters$Content)
k <- stringdot( "spectrum", length = n ) 
K <- kernelMatrix(k, reuters$Content)
R <- cluster(K, clusters, TRUE)
writeInFile("SPECTRUM KERNEL", R)
suppressWarnings(write( as.character(n),file = 'ResultsClustering.txt', append=TRUE))
##Calculate the clusteres using the constant kernel
k <- stringdot("constant",length = 2)
K <- kernelMatrix(k, reuters$Content)
R <- cluster(K, clusters, TRUE)
writeInFile("CONSTANT KERNEL",R)
suppressWarnings(write( as.character(n),file = 'ResultsClustering.txt', append=TRUE))
##Calculate the clusters using the bound range kernel
n = bestLength("boundrange", clusters, reuters$Content)
k <- stringdot("boundrange",length=n)
K <- kernelMatrix(k, reuters$Content)
R <- cluster(K, clusters, TRUE)
writeInFile("BOUNDRANGE KERNEL",R)
suppressWarnings(write( as.character(n),file = 'ResultsClustering.txt', append=TRUE))
##Calculate the cluster using the exponential kernel
##Define lambda range
lambdas = seq(1.2,2,0.2)
optimizeExponential <-function(data, k, lambdas){
  max <- 10#maxLength(data)
  Scores <- matrix(NA, length(lambdas), max-1)
  for(i in 2:max){
    for(j in 1:length(lambdas)){
      param <- stringdot(length = i, lambda = lambdas[j], type="exponential" )
      matrix <- kernelMatrix( param, data)
      possibleError <- tryCatch(
        sc <- kkmeans(matrix, k, nstart = 10),
        error=function(e) print("Param combination does not work")
      )
      if(!inherits(possibleError, "error")){
        Scores[j,i-1] = sum(withinss(sc))
      }
    }
  }
  which(Scores == min(Scores, na.rm = TRUE), arr.ind = TRUE) 
}
m <- optimizeExponential(reuters$Content, clusters, lambdas )
k <- stringdot(length = m[2]+1, lambda = lambdas[m[1]],type = "exponential" )
K <- kernelMatrix(k, reuters$Content)
R<-  cluster(K, clusters, TRUE)
writeInFile("EXPONENTIAL KERNEL",R)
suppressWarnings(write( paste(as.character(m[2]+1),as.character(lambdas[m[1]])," "),file = 'ResultsClustering.txt', append=TRUE))
##for gap Kernel
optimizeGAP <- function(data, k, lambdas, begin, limit){
  print("enter gap")
  Scores <- matrix(NA, length(lambdas), limit)
  for(i in begin:limit){
    for(j in 1:length(lambdas)){
      param <- makeGapKernel( lambdas[j], i) 
      matrix <- kernelMatrix( param, data)
      possibleError <- tryCatch({
        sc <- kkmeans(matrix, k, nstart = 10)
        Scores[j,i-1] = sum(withinss(sc))},
        error=function(e) print("param combination does not work")
      )
    }
  }
  which(Scores == min(Scores, na.rm = TRUE), arr.ind = TRUE) 
}
##Launch the gap kernel
lambdas = seq(0.1,0.9,0.1)
fc = optimizeGAP(reuters$Content,clusters,lambdas,6,8)
k <- makeGapKernel(lambdas[fc[1]], fc[2]+1)
K <- kernelMatrix(k, reuters$Content)
writeInFile("GAP KERNEL",R)
suppressWarnings(write( paste(as.character(fc[2]+1),as.character(lambdas[fc[1]])," "),file = 'ResultsClustering.txt', append=TRUE))
##Clustering using the connect kernel
k <- new("kernel", .Data=connect, kpar=list())
K <- kernelMatrix(k ,reuters$Content)
R <- cluster(K, clusters, TRUE)
writeInFile("CONNECT",R)
