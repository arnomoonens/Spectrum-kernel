source("coSequenceCPP.R")
source("gapkernel.R")

data <- read.csv("data.csv", stringsAsFactors = FALSE)
data <- data[-1]  # Remove first column (indices)
data <- data[sample(nrow(data), round(0.05 * nrow(data))),]
data <- subset(data, target<2)  # Only use three different target values

data$target <- factor(data$target)

# Total number of training examples
length(data$target)

# Number of training examples per topic
table(data$target)

(N <- dim(data)[1])

# Shuffle the data
set.seed(6)
data <- data[sample(1:N, N),]

#Calculate the Kernels
myKernel <- new("kernel",.Data=coSequenceKernelCPP,kpar=list())
start.time <- Sys.time()
K2 <- kernelMatrix(myKernel,data$text)
end.time <- Sys.time()
time.taken <- end.time - start.time
t<-time.taken

k <- stringdot("spectral",length=6,normalized = TRUE)
K <- kernelMatrix(k, data$text)
k3 <- makeCppKernel(0.7, 5)
start.time <-Sys.time()
K3 <- kernelMatrix(k3, data$text)
end.time <- Sys.time()
time.taken <- end.time - start.time
t<-time.taken
ntrain <- round(N*2/3)     # number of training examples
tindex <- sample(N,ntrain) # indices of training examples

## The fit a SVM in the train part
svm.train <- ksvm (K[tindex,tindex],data$target[tindex], type="C-svc", kernel='matrix')
svm.train2 <- ksvm (K2[tindex,tindex],data$target[tindex], type="C-svc", kernel='matrix')
svm.train3 <- ksvm (K3[tindex,tindex],data$target[tindex], type="C-svc", kernel='matrix')

# First the test-vs-train matrix
testK <- K[-tindex,tindex]
testK2 <-K[-tindex,tindex]
testK3 <-K[-tindex,tindex]
# then we extract the SV from the train
testK <- testK[,SVindex(svm.train),drop=FALSE]
testK2 <- testK2[,SVindex(svm.train2),drop=FALSE]
testK3 <- testK3[,SVindex(svm.train3),drop=FALSE]
# Now we can predict the test data
# Warning: here we MUST convert the matrix testK to a 'kernelMatrix'
y1 <- predict(svm.train,as.kernelMatrix(testK))
y2 <- predict(svm.train2,as.kernelMatrix(testK2))
y3 <- predict(svm.train3,as.kernelMatrix(testK3))

table (pred=y1, truth=data$target[-tindex])
table (pred=y2, truth=data$target[-tindex])
table (pred=y3, truth=data$target[-tindex])

cat('Error rate = ',100*sum(y1!=data$target[-tindex])/length(y1),'%')
cat('Error rate = ',100*sum(y2!=data$target[-tindex])/length(y2),'%')
cat('Error rate = ',100*sum(y3!=data$target[-tindex])/length(y3),'%')





