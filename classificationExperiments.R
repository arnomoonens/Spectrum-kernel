library(kernlab)
library(tm)
source("gapkernel.R")

set.seed(6)
data <- read.csv("data.csv", stringsAsFactors = FALSE)
data <- data[-1]  # Remove first column (indices)
data <- subset(data, target<2)  # Only use three different target values
data <- data[sample(nrow(data), round(0.15 * nrow(data))),]

docs <- Corpus(VectorSource(data$text))
docs <- tm_map(docs, removePunctuation)
docs <- tm_map(docs, stripWhitespace)
docs <- tm_map(docs, content_transformer(tolower))
data$text <- sapply(docs, function(x){x$content})

data$text <- substr(data$text, 1, 750)  # Limit number of characters per text

data$target <- factor(data$target)

(N <- nrow(data))
train.indices <- sample(N,round(4*N/5))

# Total number of data points
length(data$target)

# Number of training examples per topic
table(data$target)

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

data.train <- data[train.indices,]
# Cross-validation for subsequence length
k <- 10 # 10-fold cross-validation
folds.indices <- 1:k
N.train <- NROW(data.train)
s <- sample(rep(1:k, length=N.train), N.train, replace=FALSE) # Say for each data point in which fold it needs to go

subseq.lengths <- c(2:6)
lambdas <- c(0.3, 0.5, 0.7, 1.0)
configs <- matrix(nrow=0, ncol=2, dimnames=list(c(), c("lambda", "p")))
mean.errors <- c()
for (subseq.length in subseq.lengths) {
  for(lambda in lambdas) {
    configs <- rbind(configs, c(lambda, subseq.length))
    kernel <- makeCppKernel(lambda, subseq.length)
    K <- kernelMatrix(kernel, data.train$text)
    cv.errors <- c()
    for (idx in folds.indices) {
      train.K <- K[s != idx,s != idx]
      train.y <- data.train[s != idx,]$target
      
      model <- ksvm(train.K, train.y, type="C-svc", scaled=c(), kernel="matrix")
      
      test.K <- K[s == idx,s != idx]
      test.K <- test.K[, SVindex(model), drop=FALSE]
      test.y <- data.train[s == idx,]$target
      
      preds <- predict(model, as.kernelMatrix(test.K))
      error.cv <- calculate.error(preds, test.y)
      cv.errors <- c(cv.errors, error.cv)
    }
    m <- mean(cv.errors)
    mean.errors <- c(mean.errors, m)
    cat("Combination done")
  }
}
min.error.idx <- which.min(mean.errors)
best.config <- configs[min.error.idx,]
best.subseq.length <- best.config[["p"]]
best.lambda <- best.config[["lambda"]]
cat("Best subsequence length =", best.subseq.length)
cat("Best lambda =", best.lambda)

data.test <- data[-train.indices,]

kernel = makeCppKernel(best.lambda, best.subseq.length)
model <- ksvm(data.train$text, data.train$target, type="C-svc", scaled=c(), kernel=kernel)

predictions <- predict(model, data.test$text)
error.test <- calculate.error(predictions, data.test$target)
cat("Error on test set:", error.test)
