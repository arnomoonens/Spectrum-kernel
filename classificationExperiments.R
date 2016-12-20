library(kernlab)
library(tm)
source("gapkernel.R")
source("coSequenceCPP.R")

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
data.test <- data[-train.indices,]

# Cross-validation for subsequence length
k <- 10 # 10-fold cross-validation
folds.indices <- 1:k
N.train <- NROW(data.train)
s <- sample(rep(1:k, length=N.train), N.train, replace=FALSE) # Say for each data point in which fold it needs to go

# Gap kernel
subseq.lengths <- c(2:6)
lambdas <- c(0.3, 0.5, 0.7, 1.0)
configs <- matrix(nrow=0, ncol=2, dimnames=list(c(), c("lambda", "p")))
gap.mean.errors <- c()
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
    gap.mean.errors <- c(gap.mean.errors, m)
    cat("Combination done")
  }
}
min.error.idx <- which.min(gap.mean.errors)
gap.best.config <- gap.configs[min.error.idx,]
gap.best.subseq.length <- gap.best.config[["p"]]
gap.best.lambda <- gap.best.config[["lambda"]]
cat("Best subsequence length =", gap.best.subseq.length)
cat("Best lambda =", gap.best.lambda)

gap.kernel = makeCppKernel(gap.best.lambda, gap.best.subseq.length)
gap.model <- ksvm(data.train$text, data.train$target, type="C-svc", scaled=c(), kernel=gap.kernel)

gap.predictions <- predict(gap.model, data.test$text)
gap.error.test <- calculate.error(gap.predictions, data.test$target)
cat("Error on test set:", gap.error.test)

# Spectrum kernel
substring.lengths <- c(2:6)
spectrum.mean.errors <- c()
for (substring.length in substring.lengths) {
  kernel <- stringdot("spectrum", length = substring.length)
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
  spectrum.mean.errors <- c(spectrum.mean.errors, m)
}
min.error.idx <- which.min(spectrum.mean.errors)
spectrum.best.substring.length <- substring.lengths[[min.error.idx]]
cat("Best substring length =", spectrum.best.substring.length)

spectrum.kernel = stringdot("spectrum", length = spectrum.best.substring.length)
spectrum.model <- ksvm(data.train$text, data.train$target, type="C-svc", scaled=c(), kernel=spectrum.kernel)

spectrum.predictions <- predict(spectrum.model, data.test$text)
spectrum.error.test <- calculate.error(spectrum.predictions, data.test$target)
cat("Error on test set:", spectrum.error.test)

# Boundrange kernel
substring.lengths <- c(2:6)
boundrange.mean.errors <- c()
for (substring.length in substring.lengths) {
  kernel <- stringdot("boundrange", length = substring.length)
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
  boundrange.mean.errors <- c(boundrange.mean.errors, m)
}
min.error.idx <- which.min(boundrange.mean.errors)
boundrange.best.substring.length <- substring.lengths[[min.error.idx]]
cat("Best substring length =", boundrange.best.substring.length)

boundrange.kernel = stringdot("boundrange", length = boundrange.best.substring.length)
boundrange.model <- ksvm(data.train$text, data.train$target, type="C-svc", scaled=c(), kernel=boundrange.kernel)

boundrange.predictions <- predict(boundrange.model, data.test$text)
boundrange.error.test <- calculate.error(boundrange.predictions, data.test$target)
cat("Error on test set:", boundrange.error.test)

# Constant kernel
substring.lengths <- c(2:6)
const.mean.errors <- c()
for (substring.length in substring.lengths) {
  kernel <- stringdot("constant", length = substring.length)
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
  const.mean.errors <- c(const.mean.errors, m)
}
min.error.idx <- which.min(const.mean.errors)
const.best.substring.length <- substring.lengths[[min.error.idx]]
cat("Best substring length =", const.best.substring.length)

const.kernel = stringdot("constant", length = const.best.substring.length)
const.model <- ksvm(data.train$text, data.train$target, type="C-svc", scaled=c(), kernel=const.kernel)

const.predictions <- predict(const.model, data.test$text)
const.error.test <- calculate.error(const.predictions, data.test$target)
cat("Error on test set:", const.error.test)

# Exponential kernel
lambdas <- c(1.1, 1.3, 1.5, 1.7, 1.9)
exp.mean.errors <- c()
for (lambda in lambdas) {
  kernel <- stringdot("exponential", lambda = lambda, length = 2)
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
  exp.mean.errors <- c(exp.mean.errors, m)
}
min.error.idx <- which.min(exp.mean.errors)
exp.best.lambda <- lambdas[[min.error.idx]]
cat("Best lambda =", exp.best.lambda)

exp.kernel = stringdot("exponential", lambda = exp.best.lambda)
exp.model <- ksvm(data.train$text, data.train$target, type="C-svc", scaled=c(), kernel=exp.kernel)

exp.predictions <- predict(exp.model, data.test$text)
exp.error.test <- calculate.error(exp.predictions, data.test$target)
cat("Error on test set:", exp.error.test)

# Cosequence kernel
coseq.kernel <- new("kernel", .Data=coSequenceKernelCPP, kpar=list())
coseq.model <- ksvm(data.train$text, data.train$target, type="C-svc", scaled=c(), kernel=coseq.kernel)

coseq.predictions <- predict(coseq.model, data.test$text)
coseq.error.test <- calculate.error(coseq.predictions, data.test$target)
cat("Error on test set:", coseq.error.test)
