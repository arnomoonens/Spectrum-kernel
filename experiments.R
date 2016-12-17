library(kernlab)
library(tm)
source("gapkernel.R")

data <- read.csv("data.csv", stringsAsFactors = FALSE)
data <- data[-1]  # Remove first column (indices)
data <- data[sample(nrow(data), round(0.1 * nrow(data))),]
data <- subset(data, target<2)  # Only use three different target values

docs <- Corpus(VectorSource(data$text))
docs <- tm_map(docs, removePunctuation)
docs <- tm_map(docs, stripWhitespace)
data$text <- sapply(docs, function(x){x$content})

data$text <- substr(data$text, 1, 750)  # Limit number of characters per text

data$target <- factor(data$target)

# Total number of training examples
length(data$target)

# Number of training examples per topic
table(data$target)

(N <- dim(data)[1])

# Shuffle the data
set.seed(6)
data <- data[sample(1:N, N),]


#k <- stringdot("spectrum", length=5, normalized=T)
# k <- makeCppKernel(0.5, 5)
# K <- kernelMatrix(k, data$text)
# 
# ntrain <- round(N*2/3)     # number of training examples
# tindex <- sample(N,ntrain) # indices of training examples
# 
# ## The fit a SVM in the train part
# svm.train <- ksvm (K[tindex,tindex],data$target[tindex], type="C-svc", kernel='matrix')
# 
# # First the test-vs-train matrix
# testK <- K[-tindex,tindex]
# # then we extract the SV from the train
# testK <- testK[,SVindex(svm.train),drop=FALSE]
# 
# # Now we can predict the test data
# # Warning: here we MUST convert the matrix testK to a 'kernelMatrix'
# y1 <- predict(svm.train,as.kernelMatrix(testK))
# 
# table (pred=y1, truth=data$target[-tindex])
# 
# cat('Error rate = ',100*sum(y1!=data$target[-tindex])/length(y1),'%')

calculate.error  <- function(predictions, true.values) {sum(predictions != true.values) / NROW(true.values)}

# Cross-validation for subsequence length
k <- 10 # 10-fold cross-validation
folds.indices <- 1:k
N <- NROW(data)
s <- sample(rep(1:k, length=N), N, replace=FALSE) # Say for each data point in which fold it needs to go

subseq.lengths <- c(2:6)
mean.errors <- c()
for (subseq.length in subseq.lengths) {
  kernel <- makeCppKernel(0.5, subseq.length)
  K <- kernelMatrix(kernel, data$text)
  cv.errors <- c()
  for (idx in folds.indices) {
    train.K <- K[s != idx,s != idx]
    train.y <- data[s != idx,]$target
    
    model <- ksvm(train.K, train.y, type="C-svc", scaled=c(), kernel="matrix")
    
    test.K <- K[s == idx,s != idx]
    test.K <- test.K[, SVindex(model), drop=FALSE]
    test.y <- data[s == idx,]$target
    
    preds <- predict(model, as.kernelMatrix(test.K))
    error.cv <- calculate.error(preds, test.y)
    cv.errors <- c(cv.errors, error.cv)
  }
  m <- mean(cv.errors)
  mean.errors <- c(mean.errors, m)
}
min.error.idx <- which.min(mean.errors)
best.subseq.length <- subseq.lengths[[min.error.idx]]


