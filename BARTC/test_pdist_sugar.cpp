#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector pdistC2(double x, NumericVector ys) {
  return pow((x - ys), 2);
}

