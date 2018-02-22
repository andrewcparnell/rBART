#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List fastLm(NumericVector yr, NumericMatrix Xr) {

   int n = Xr.nrow(), k = Xr.ncol();
   
   arma::mat X(Xr.begin(), n, k, false); 
   arma::colvec y(yr.begin(), yr.size(), false);
   
   arma::colvec coef = arma::solve(X, y); 
   arma::colvec resid = y - X*coef; 
   
   double sig2 = arma::as_scalar(arma::trans(resid)*resid/(n-k));
   arma::colvec stderrest = arma::sqrt(
       sig2 * arma::diagvec( arma::inv(arma::trans(X)*X)) );
   
   return List::create(Named("coefficients") = coef,
                       Named("stderr")       = stderrest);
}