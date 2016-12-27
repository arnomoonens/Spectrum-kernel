## This part asks for the discovering of the
## correct number of clusters.
library(kernlab)
library(ggplot2)
require(plot3D)
library(tm)
source("connect.R")
source("gapkernel.R")

findBestK <- function(K, limit,name, kReal){
  ss <-vector(length = limit)
  set.seed(20)
  for( i in seq(2,limit+1)) {
          possibleError <- tryCatch({
            ss[i-1] <- sum(withinss(kkmeans(K, centers = i, nstart = 10)))
            },
            error=function(e) print("Param combination does not work")
          )
  }
  k = which.min(ss)
  jpeg(paste("charts/",name,".jpg"))
  plot(seq(2,limit+1), ss,  type="o",
       main = paste(name,"-kernel ", "k real = ",as.character(kReal) ),
       xlab = "different k", ylab="sum of ss")
  dev.off()
}
##Size of the samples
N = 30
## Read table
reuters <- read.table("reuters.txt.gz", header=TRUE)
reuters$Content <- as.character(reuters$Content)  
reuters$Topic <- factor(reuters$Topic)
## Sample a number cat of categories
cat = 15
set.seed(6)
topics <- factor( sample(unique(reuters$Topic),cat) )
reuters <- reuters[ sapply(reuters[,1], function (x) x %in% topics) , ]
set.seed(6)
reuters <- reuters[sample(1:nrow(reuters),N),]
##Preprocess the data
docs <- Corpus(VectorSource(reuters$Content))
docs <- tm_map(docs, removePunctuation)
docs <- tm_map(docs, stripWhitespace)
docs <- tm_map(docs,removeNumbers)
reuters$Content <- sapply(docs, function(x){x$content})
##Correct number of categories
kReal = length(unique(reuters[,1]))


limit = 11

k <- stringdot("spectrum", length=2)
K <- kernelMatrix(k, reuters$Content)
findBestK(K,limit,"SPECTRUM",kReal ) 

k <- new("kernel", .Data=connect, kpar=list())
K <- kernelMatrix(k ,reuters$Content)
findBestK(K,limit,"CONNECT",kReal ) 

k <- makeGapKernel(0.1, 8)
K <- kernelMatrix(k, reuters$Content)
findBestK(K,limit,"GAP",kReal ) 

k <- stringdot("exponential",length= 4, lambda = 1.4 )
K <- kernelMatrix(k, reuters$Content)
findBestK(K, limit,"EXPONENTIAL",kReal ) 



