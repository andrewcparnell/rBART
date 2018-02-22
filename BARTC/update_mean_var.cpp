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
//############################################################################

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List update_mean_var(NumericMatrix tree_table,IntegerMatrix tree_matrix,NumericVector resids,double a,double sigma,double mu_mu){
  List update_params(2);
  IntegerVector terminal_nodes;
   
  NumericVector term_obs;
  
 /* for(int i=0;i<tree_table.nrow();i++){
    if(tree_table(i,4)==-1){
      terminal_nodes.push_back(i+1);
    }
  }*/
  terminal_nodes= find_term_nodes(tree_table);
  
  NumericVector mu_ij(terminal_nodes.size());
  NumericVector Tj(terminal_nodes.size());
  NumericVector new_mean(terminal_nodes.size());
  NumericVector new_var(terminal_nodes.size());
  
  for(int k=0;k< terminal_nodes.size();k++){
    //update the node means    
    //find which observations are assigned to terminal node k in tree i
     
    term_obs=find_term_obs(tree_matrix,terminal_nodes[k]);
    
    //get the number of observations in node k
    
    Tj[k]=term_obs.size();
      
    for(int i=0;i<term_obs.size();i++){
      mu_ij[k]+=resids[term_obs[i]];
    } 
      
    new_mean[k]=(mu_ij[k]+a*mu_mu)/(Tj[k]+a);;
  
    new_var[k]=(1/((1/sigma)*(Tj[k]+a)));
    
    NumericVector temp;
    term_obs=temp;

  }  
  
  //only thing changing is mu_ij the terminal node mean
  //everything else was set outside the loop
  
  update_params[0]=new_mean;
  update_params[1]=new_var;
  
  return(update_params);
}