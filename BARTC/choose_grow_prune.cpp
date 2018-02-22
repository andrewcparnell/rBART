#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

double get_tree_prior(NumericMatrix tree_table,IntegerMatrix tree_matrix,double alpha,double beta){
  double propsplit=1;
  IntegerVector d;
  IntegerVector int_nodes_index;
  NumericVector internal_nodes_prop;
  
  for(int i=0;i<tree_table.nrow();i++){
    if(tree_table(i,4)==1){
      internal_nodes_prop.push_back(i+1); 
    }
  }
  
  for(int k=0;k<internal_nodes_prop.size();k++){ 
    for(int i=0;i<tree_matrix.nrow();i++){
      for(int j=0;j<tree_matrix.ncol();j++){
        
        if(tree_matrix(i,j)==internal_nodes_prop[k]){
          int_nodes_index.push_back(j+1);  
        }
      }
    }

    if(int_nodes_index.size()!=0){
      d=unique(int_nodes_index);
      double d1=d[0];
      propsplit*=alpha*pow((d1+1),-beta) ;
      
    }
    IntegerVector temp;
    int_nodes_index=temp;
  } 
  return(propsplit);
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

double choose_grow_prune(NumericMatrix proposal_table,IntegerMatrix proposal_matrix,NumericMatrix prior_table,IntegerMatrix prior_matrix,double alpha,double beta){
   double tree_ratio;
   double propsplit;
   double prisplit;
   
   propsplit=get_tree_prior(proposal_table,proposal_matrix,alpha,beta);
  
  //and do the same for the prior tree
  prisplit=get_tree_prior(prior_table,prior_matrix,alpha,beta);
  tree_ratio=propsplit/prisplit;
  return(tree_ratio);
}