#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

NumericVector find_term_obs(IntegerMatrix tree_matrix,int terminal_node){
  NumericVector term_obs;
for(int i=0;i<tree_matrix.nrow();i++){
      for(int j=0;j<tree_matrix.ncol();j++){
        
        if(tree_matrix(i,j)==terminal_node){
          term_obs.push_back(i);  
        }
      }
    }
    return(term_obs);
}