library(kernlab)
set.seed(6)
n = 30
reuters <- read.table("reuters.txt.gz", header=T)
# We leave only three topics for analysis: Crude Oil, Coffee and Grain-related news
reuters <- reuters[reuters$Topic == "crude" | reuters$Topic == "grain" | reuters$Topic == "coffee",]
reuters <- reuters[sample(1:nrow(reuters),n),]
reuters$Content <- as.character(reuters$Content)    # R originally loads this as factor, so needs fixing
reuters$Topic <- factor(reuters$Topic)              # re-level the factor to have only three levels
levels(reuters$Topic)
length(reuters$Topic)
table(reuters$Topic)

getClusters <- function(kernelMatrix){
  sc <- kkmeans(K, reuters$Topic,centers=3, nstart = 20)
  print( table(sc, reuters$Topic) )
  print ( withinss(sc) )
}
#Calculate the Kernels
#myKernel <- new("kernel",.Data=coSequenceKernelCPP,kpar=list())
#K2 <- kernelMatrix(myKernel,data$text)
#k3 <- makeCppKernel(0.7, 5)
set.seed(20)
k <- stringdot("spectral",length=6)
K <- kernelMatrix(k, reuters$Content)
getClusters(K)
k <- stringdot("boundrange",length=6)
K <- kernelMatrix(k, reuters$Content)
getClusters(K)
k <- stringdot("constant",length=6)
K <- kernelMatrix(k, reuters$Content)
getClusters(K)
k <- stringdot("exponential",length=6,lambda = 1.1)
K <- kernelMatrix(k, reuters$Content)
getClusters(K)





