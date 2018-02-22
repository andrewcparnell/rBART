#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
IntegerMatrix set_daughter(int left_daughter,int right_daughter,IntegerVector ld_obs,IntegerVector rd_obs,IntegerMatrix tree_matrix,double term_cols){
  for(int i=0;i<ld_obs.size();i++){
    tree_matrix(ld_obs[i],term_cols+1)=left_daughter;    
  }
  for(int i=0;i<rd_obs.size();i++){
    tree_matrix(rd_obs[i],term_cols+1)=right_daughter;    
  }    
  return(tree_matrix);
}