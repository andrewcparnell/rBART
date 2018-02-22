#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix start_predy(int num_trees,int n,double start_mean,double start_sd){
  NumericMatrix predy(n,num_trees);
  int ncol=predy.ncol();
  int nrow=predy.nrow();
  for(int i=0;i<ncol;i++){
    for(int j=0;j<nrow;j++){
      predy(j,i)=R::rnorm(start_mean,start_sd);
    }  
  }
  return(predy);
}