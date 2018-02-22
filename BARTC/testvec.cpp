#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector testvec(){
  double rand=R::rnorm(start_mean,start_sd);
NumericVector testrow = NumericVector::create(0,0,0,0,-1,rand,0);
return testrow;
}