#include <Rcpp.h>
#include <math.h> // pow
using namespace Rcpp;

// [[Rcpp::export]]
long double cppKernel(std::string x, std::string y, double lambda, unsigned int p) {
  int r = x.length();
  int c = y.length();
  unsigned int i = 0, j = 0;
  long double **A = new long double* [r];
  long double **B = new long double* [r+1];
  
  B[0] = new long double[c+1];
  for(i = 0; i < r; i++) {
    A[i] = new long double[c];
    B[i+1] = new long double[c+1];
    for(j = 0; j < c; j++) {
      B[i][j] = 0;
      B[i+1][j+1] = 0;
      A[i][j] = (x[i] == y[j] ? pow(lambda, 2) : 0);
    }
  }
  long double* k = new long double[p];
  k[0] = 0;
  for(int l = 1; l < p; l++) {
    k[l] = 0;
    for(i = 0; i < r; i++) {
      for(j = 0; j < c; j++) {
        B[i+1][j+1] = A[i][j] + lambda * B[i][j+1] + lambda * B[i+1][j] - pow(lambda, 2) * B[i][j];
        if(x[i] == y[j]) {
          A[i][j] = pow(lambda, 2) * B[i][j];
          k[l] = k[l] + A[i][j];
        }
      }
    }
  }
  delete[] B[0];
  for(i = 0; i < r; i++) {
    delete[] A[i];
    delete[] B[i+1];
  }
  delete[] A;
  delete[] B;
  long double result = k[p-1];
  delete[] k;
  return result;
}


// (test) code to run after compilation
/*** R
x <- "dfqsdf qdsf qsdfqds qsdf  dfqsdfdf qdsf qsqsd f qsd f qsd fd qsqsd f qsd f qsd fd df qdsf qsqsd f qsd f qsd fd qsdf qsdf qdsf qsqsd f qsd f qsd fd qsdf  dfqsdf qdsf qsdfqds f qsd f qsd f qsd fd qsdf  dfqsdf qdsf qsdfqds f qsd f qsd f qsd fd qsdf  dfqsdf qdsf qsdfqds f qsd f qsd f qsd fd qsdf  dfqsdf qdsf qsdfqds f qsd f qsd f qsd fd qsdf  dfqsdf qdsf qsdfqds f qsd f qsd f qsd fd qsdf  dfqsdf qdsf qsdfqds f qsd f qsd f qsd fd qsdf  dfqsdf qdsf qsdfqds f qsd f qsd f qsd fd qsdf  dfqsdf qdsf qsdfqds f qsd f qsd f qsd fd qsdf  dfqsdf qdsf qsdfqds f qsd f qsd f qsd fd qsdf  "
cppKernel(x, x, 0.7, 2)
*/
