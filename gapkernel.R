# Gap weighted subsequence kernel
# By Arno Moonens

gapkernel <- function(lambda, p) {
  gapkernelfunction <- function(x, y) 
  {
    r <- nchar(x)
    c <- nchar(y)
    A <- matrix(nrow = r+1, ncol = c+1)
    B <- matrix(nrow = r+1, ncol = c+1)
    for(i in c(1:r)) {
      B[i+1, 1] <- 0
      for(j in c(1:c)) {
        A[i+1, j+1] <- 0
        if(substr(x, i, i) == substr(y, j, j)) {
          A[i+1, j+1] <- lambda^2
        }
      }
    }
    B[1,1] <- 0
    for(j in c(1:c)) {
      B[1, j+1] <- 0
    }
    k <- c(1:p)
    for(l in c(2:p)) {
      k[l] <- 0
      for(i in c(1:(r-1))) {
        for(j in c(1:(c-1))) {
          B[i+1, j+1] <- A[i+1, j+1] + lambda * B[i, j+1] + lambda * B[i+1, j] + lambda^2 * B[i, j]
          if(substr(x, i, i) == substr(y, j, j)) {
            A[i + 1, j + 1] <- lambda ^ 2 * B[i, j]
            k[l] <- k[l] + A[i + 1, j + 1]
          }
        }
      }
    }
    return(k[p])
  }
  return(new("kernel",.Data=gapkernelfunction,kpar=list(lambda = lambda, p = p)))
}

# Faster version
gapkernel2 <- function(lambda, p) {
  gapkernelfunction <- function(x, y) 
  {
    x <- strsplit(x, split="")[[1]]
    y <- strsplit(y, split="")[[1]]
    cat("splitting done")
    r <- length(x)
    c <- length(y)
    A <- matrix(0, nrow = r+1, ncol = c+1)
    B <- matrix(0, nrow = r+1, ncol = c+1)
    for(i in c(1:r)) {
      for(j in c(1:c)) {
        if(x[i] == y[j]) {
          A[i+1, j+1] <- lambda^2
        }
      }
    }
    cat("First loop done")
    k <- c(1:p)
    for(l in c(2:p)) {
      k[l] <- 0
      for(i in c(1:(r-1))) {
        for(j in c(1:(c-1))) {
          B[i+1, j+1] <- A[i+1, j+1] + lambda * B[i, j+1] + lambda * B[i+1, j] + lambda^2 * B[i, j]
          if(x[i] == y[j]) {
            A[i + 1, j + 1] <- lambda ^ 2 * B[i, j]
            k[l] <- k[l] + A[i + 1, j + 1]
          }
        }
      }
    }
    return(k[p])
  }
  return(new("kernel",.Data=gapkernelfunction,kpar=list(lambda = lambda, p = p)))
}

Rcpp::sourceCpp('gapkernel.cpp')
makeGapKernel <- function(lambda, p) {
  f <- function(x, y) {
    return(cppKernel(x, y, lambda, p))
  }
  return(new("kernel",.Data=f,kpar=list(lambda = lambda, p = p)))
}
