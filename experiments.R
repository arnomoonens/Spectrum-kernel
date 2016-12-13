library(kernlab)
source("gapkernel.R")

data <- read.csv("data.csv", stringsAsFactors = FALSE)
data <- data[-1]  # Remove first column (indices)
data <- data[sample(nrow(data), round(0.1 * nrow(data))),]
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


#k <- stringdot("spectrum", length=5, normalized=T)
k <- makeCppKernel(0.7, 5)
K <- kernelMatrix(k, data$text)

ntrain <- round(N*2/3)     # number of training examples
tindex <- sample(N,ntrain) # indices of training examples

## The fit a SVM in the train part
svm.train <- ksvm (K[tindex,tindex],data$target[tindex], type="C-svc", kernel='matrix')

# First the test-vs-train matrix
testK <- K[-tindex,tindex]
# then we extract the SV from the train
testK <- testK[,SVindex(svm.train),drop=FALSE]

# Now we can predict the test data
# Warning: here we MUST convert the matrix testK to a 'kernelMatrix'
y1 <- predict(svm.train,as.kernelMatrix(testK))

table (pred=y1, truth=data$target[-tindex])

cat('Error rate = ',100*sum(y1!=data$target[-tindex])/length(y1),'%')
