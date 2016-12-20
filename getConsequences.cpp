#include <Rcpp.h>
#include <math.h>       /* pow */
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericVector matches(NumericVector  x, NumericVector  y){
  int i = 0;
  int j = 0;
  NumericVector indexes;
  while(true){
    bool found = false;
    bool moved = false;
    if(x[i] + 1 == y[j]){
      indexes.push_back(y[j]);
      i++; j++;
      found = true;
      moved= true;
    }
    if(i == x.size() && j == y.size()) break;
    if(!found){
      if(x[i] < y[j] && i < x.size()) {i++; moved= true;}
      if(x[i] >= y[j] && j < y.size()) {j++; moved = true;}
    }
    if (moved ==false) break;
  }
  return indexes;
}

// [[Rcpp::export]]
List commonSubSequences(NumericVector x_1, NumericVector  x_2, NumericVector 
                          y_1, NumericVector  y_2 ){
  NumericVector matchesX1 = matches(x_1,y_1);
  NumericVector matchesX2 = matches(x_2,y_2);
  if (matchesX1.size() == 0 || matchesX2.size() == 0){
    return NULL;
  }
  return List::create(matchesX1,matchesX2);
}

// [[Rcpp::export]]
List sequences(List extendables, List index) {
  NumericVector vec = as<NumericVector>(as<List>(extendables[0])[0]);
  List k_extendables = NULL;
  int counter = 0;
  NumericVector x_1, x_2, y_1, y_2;
  List tmp;
  for(int i = 0; i < extendables.size(); i ++){
    for(int j = 0; j < index.size(); j ++){
      x_1 = as<NumericVector>(as<List>(extendables[i])[0]);
      x_2 = as<NumericVector>(as<List>(extendables[i])[1]);
      y_1 = as<NumericVector>(as<List>(index[j])[0]);
      y_2 = as<NumericVector>(as<List>(index[j])[1]);
      tmp = commonSubSequences(x_1,x_2,y_1,y_2);
      if (!tmp.isNULL() && tmp.size() != 0) {
        k_extendables.push_back(tmp);
      }
    }
  }
  return k_extendables;
}
