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