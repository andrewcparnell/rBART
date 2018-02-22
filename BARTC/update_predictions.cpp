//############################################################################
//need to move find_term_nodes to a header file
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

NumericVector find_term_nodes(NumericMatrix tree_table){
  
  NumericVector terminal_nodes;
  
  for(int i=0;i<tree_table.nrow();i++){
    if(tree_table(i,4)==-1){
      terminal_nodes.push_back(i+1);
    }
  }
  return(terminal_nodes);
}  

//############################################################################
//need to move find_term_nodes to a header file
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

//########################################################################################//
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List update_predictions(NumericMatrix tree_table,IntegerMatrix tree_matrix,NumericVector new_mean,NumericVector new_var,int n){
  List updated_preds(2);
  NumericVector new_preds(n);
  NumericVector terminal_nodes;
  NumericVector term_obs;
  
  terminal_nodes=find_term_nodes(tree_table);
  
   for(int k=0;k<terminal_nodes.size();k++){
    
    term_obs=find_term_obs(tree_matrix,terminal_nodes[k]);
            
    //update the mean of the selected tree nodes:
    tree_table(terminal_nodes[k]-1,5)= new_mean[k];
    tree_table(terminal_nodes[k]-1,6)=sqrt(new_var[k]);
    
    //update residuals for next iteration
   
    for(int i=0;i<term_obs.size();i++){
      new_preds[term_obs[i]]=new_mean[k];
    }
  }
  
  updated_preds[0]=tree_table;
  updated_preds[1]=new_preds;
  
  return(updated_preds);
}