#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector calc_rowsums(NumericMatrix predictions,NumericVector response){
  NumericVector row_sums=response.size();
  NumericVector resids=response.size();
  int ncol=predictions.ncol();
  int nrow=predictions.nrow();
  for(int j=0;j<nrow;j++){
    for(int i=0;i<ncol;i++){  
      row_sums[j]+=predictions(j,i);
    }  
  }
  return(row_sums);
}
// [[Rcpp::export]]
NumericVector calculate_resids(NumericMatrix predictions,NumericVector response){
  NumericVector resids=response.size();
  NumericVector row_sums=calc_rowsums(predictions,response);
  resids=response - row_sums;
  return(resids);
}
