#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
IntegerVector get_imp_vars(IntegerVector tree_table,int num_col,IntegerVector current_vars){
  
  IntegerVector vars_chosen=sort_unique(tree_table);
  
  if(vars_chosen[0]==0){
    vars_chosen.erase(0);
  }
  
  for(int i=0;i<tree_table.size();i++){
    
    if(tree_table[i]!=0){      
      current_vars[tree_table[i]-1]+=1;
    }
    
  }
  
  return(current_vars);
}