#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector testvec(){
NumericVector testrow = NumericVector::create(0,0,0,0,-1,2,0);
return testrow;
}