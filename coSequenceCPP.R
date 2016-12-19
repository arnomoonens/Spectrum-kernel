library(kernlab)
library(Matrix)

Rcpp::sourceCpp('getConsequences.cpp')
kernel<- function(x,y){
  index.x <- vector(mode="list", length = 27)
  index.y <- vector(mode="list", length = 27)
  get <- function(x,i) {substr(x,i,i)}
  map <- function(x) { match( x, letters, nomatch = 27 ) }#use sapply
  #mapply(function(x,y) {index.x[[map(x)]] = c(index.x[[ map(x) ]], y ) }, x, 1:nchar(x) )
  #mapply(function(x,y) {index.x[[map(x)]] = c(index.x[[ map(x) ]], y ) }, y, 1:nchar(y) )
  
  for(i in 1:nchar(x)) { index.x[[map(get(x,i))]] = c(index.x[[map(get(x,i))]],i) }
  for(i in 1:nchar(y)) { index.y[[map(get(y,i))]] = c(index.y[[map(get(y,i))]],i) }
  index.xy = list()
  j <-1
  for(i in 1:26){
    if(!is.null(index.x[[i]]) & !is.null(index.y[[i]])){
      index.xy[[j]] = list(index.x[[i]],index.y[[i]])
      j <- j+1
    }
  }
  extendables = index.xy
  NumberCoSequences = list(vector(mode="list", length=min(nchar(x),nchar(y))))
  for(k in 1:min(nchar(x),nchar(y))){
    NumberCoSequences[[k]] = sum( sapply(extendables, function(x) { length(x[[1]])* length(x[[2]]) } ) )
    extendables <- sequences(extendables, index.xy)
    if( length(extendables)==0) break
  }
  #NumberCoSequences
  as.double(sum(mapply(sum,NumberCoSequences)))
}

coSequenceKernelCPP <- function(x,y){
  as.double( kernel(x,y) ) / sqrt( kernel(x,x)* kernel(y,y) )
}
