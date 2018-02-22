#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List start_matrix(int num_trees,int n){
  List prior_matrix(num_trees);
  NumericMatrix mat(n,1);
  std::fill(mat.begin(), mat.end(), 1);
  for(int i=0; i<prior_matrix.size();i++){
     prior_matrix[i]=mat;
  }
  return(prior_matrix);
}