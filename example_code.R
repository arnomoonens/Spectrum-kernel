set.seed(6046)
library(kernlab)

################################################
# SVM for classification with a string kernel
################################################

## We are going to use a slightly-processed version of the famous
## Reuters news articles dataset.  All articles with no Topic
## annotations are dropped. The text of each article is converted to
## lowercase, whitespace is normalized to single-spaces.  Only the
## first term from the Topic annotation list is retained (some
## articles have several topics assigned).  

## The resulting dataset is a list of pairs (Topic, News content) We willl only use three topics for analysis: Crude Oil, Coffee and Grain-related news

## The resulting data frame contains 994 news items on crude oil,
## coffee and grain. The news text is the column "Content" and its
## category is the column "Topic". The goal is to create a classifier
## for the news articles.

## Note that we can directly read the compressed version (reuters.txt.gz). 
## There is no need to unpack the gz file; for local files R handles unpacking automagically

reuters <- read.table("reuters.txt.gz", header=T)

# We leave only three topics for analysis: Crude Oil, Coffee and Grain-related news
reuters <- reuters[reuters$Topic == "crude" | reuters$Topic == "grain" | reuters$Topic == "coffee",]

reuters$Content <- as.character(reuters$Content)    # R originally loads this as factor, so needs fixing
reuters$Topic <- factor(reuters$Topic)              # re-level the factor to have only three levels

levels(reuters$Topic)

length(reuters$Topic)

table(reuters$Topic)

## an example of a text about coffee
reuters[2,]

## an example of a text about grain
reuters[7,]

## an example of a text about crude oil
reuters[12,]

(N <- dim(reuters)[1])  # number of rows

# we shuffle the data first
set.seed(12)
reuters <- reuters[sample(1:N, N),]

# To deal with textual data we need to use a string kernel. Several such kernels are implemented in the "stringdot" method of the kernlab package. We shall use the simplest one: the p-spectrum kernel. The feature map represents the string as a multiset of its substrings of length p

# Example, for p=2 we have

# phi("ababc") = ("ab" -> 2, "ba" -> 1, "bc" --> 1, other -> 0)

# we can define a normalized 3-spectrum kernel (p is length)
k <- stringdot("spectrum", length=3, normalized=T)

# Let's see some examples:

k("I did it my way", "I did it my way")

k("He did it his way", "I did it my way")

k("I did it my way", "She did it her way")

k("I did it my way", "Let's get our way out")

## We start by doing a kPCA (we'll see this in the next lecture)

## first we define a modified plotting function 

plotting <-function (kernelfu, kerneln)
{
  xpercent <- eig(kernelfu)[1]/sum(eig(kernelfu))*100
  ypercent <- eig(kernelfu)[2]/sum(eig(kernelfu))*100
  
  plot(rotated(kernelfu), col=as.integer(reuters$Topic),
       main=paste(paste("Kernel PCA (", kerneln, ")", format(xpercent+ypercent,digits=3)), "%"),
       xlab=paste("1st PC -", format(xpercent,digits=3), "%"),
       ylab=paste("2nd PC -", format(ypercent,digits=3), "%"))
}

## Create a kernel matrix using 'k' as kernel

k <- stringdot("spectrum", length=5, normalized=T)
K <- kernelMatrix(k, reuters$Content)
dim(K)

K[2,2]

K[2,3:10]

## Plot the result using the first 2 PCs (we can add colors for the two classes)

kpc.reuters <- kpca (K, features=2, kernel="matrix")
plotting (kpc.reuters,"5 - spectrum kernel")

## finally add a legend
legend("topleft", legend=c("crude oil", "coffee","grain"),    
       pch=c(1,1),                    # gives appropriate symbols
       col=c("red","black", "green")) # gives the correct color

## We can also train a SVM using this kernel matrix in the training set

## First we should split the data into learning (2/3) and test (1/3) parts
ntrain <- round(N*2/3)     # number of training examples
tindex <- sample(N,ntrain) # indices of training examples

## The fit a SVM in the train part
svm1.train <- ksvm (K[tindex,tindex],reuters$Topic[tindex], type="C-svc", kernel='matrix')

## and make it predict the test part

## Let's call SV the set of obtained support vectors

## Then it becomes tricky. We must compute the test-vs-SV kernel matrix
## which we do in two phases:

# First the test-vs-train matrix
testK <- K[-tindex,tindex]
# then we extract the SV from the train
testK <- testK[,SVindex(svm1.train),drop=FALSE]

# Now we can predict the test data
# Warning: here we MUST convert the matrix testK to a 'kernelMatrix'
y1 <- predict(svm1.train,as.kernelMatrix(testK))

table (pred=y1, truth=reuters$Topic[-tindex])

cat('Error rate = ',100*sum(y1!=reuters$Topic[-tindex])/length(y1),'%')

## now we define a 3D plotting function

library("rgl")
open3d()

plotting3D <-function (kernelfu, kerneln)
{
  xpercent <- eig(kernelfu)[1]/sum(eig(kernelfu))*100
  ypercent <- eig(kernelfu)[2]/sum(eig(kernelfu))*100
  zpercent <- eig(kernelfu)[3]/sum(eig(kernelfu))*100
  
  # resize window
  par3d(windowRect = c(100, 100, 612, 612))
  
  plot3d(rotated(kernelfu), 
         col  = as.integer(reuters$Topic),
         xlab = paste("1st PC -", format(xpercent,digits=3), "%"),
         ylab = paste("2nd PC -", format(ypercent,digits=3), "%"),
         zlab = paste("3rd PC -", format(zpercent,digits=3), "%"),
         main = paste("Kernel PCA"), 
         sub = "red - crude oil | black - coffee | green - grain",
         top = TRUE, aspect = FALSE, expand = 1.03)
}

kpc.reuters <- kpca (K, features=3, kernel="matrix")
plotting3D (kpc.reuters,"5 - spectrum kernel")