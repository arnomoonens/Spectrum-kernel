#include <Rcpp.h>
#include <math.h>       /* pow */
using namespace Rcpp;

// [[Rcpp::export]]
double cppKernel(const std::string x, const std::string y, double lambda, int p) {
  int r = x.length();
  int c = y.length();
  double A[r][c];
  double B[r][c];
  for(int i = 0; i < r; i++) {
    for(int j = 0; j < c; j++)
      if(x[i] == y[j]) {
        A[i][j] = pow(lambda, 2);
      }
  }
  double k[p+1];
  for(int l = 2; l <= p; l++) {
    k[l] = 0;
    for(int i = 0; i < r-1; i++) {
      for(int j = 0; j < c-1; j++) {
        B[i+1][j+1] = A[i+1][j+1] + lambda * B[i][j+1] + lambda * B[i+1][j] + pow(lambda, 2) * B[i][j];
        if(x[i] == y[j]) {
          A[i + 1][j + 1] = pow(lambda, 2) * B[i][j];
          k[l] = k[l] + A[i + 1][j + 1];
        }
      }
    }
  }
  return k[p];
}


// (test) code to run after compilation
/*** R
#cppKernel("abcde", "bcde", 0.7, 2)
*/
