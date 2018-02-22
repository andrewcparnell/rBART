#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
RcppExport SEXP getRGamma() {
RNGScope scope;
NumericVector x = rgamma( 10, 1, 1 );
return x;
}
