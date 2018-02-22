#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
  IntegerVector find_term_child(NumericMatrix tree_table){
  NumericVector node_status=tree_table.column(4);
  IntegerVector internal_nodes_prop;
  
  for(int i=0;i<node_status.size();i++){
   // std::cout<<"i is"<<" "<<i<<"node status is"<<node_status[i]<<"\n";
    if(node_status[i]==1){
    //std::cout<<"NODE STATUS=1!!"<<"\n";
      internal_nodes_prop.push_back(i+1);
      
    }
  }
  
  //std::cout<<"Internal_nodes_prop"<<internal_nodes_prop[1]<<" "<<internal_nodes_prop.size()<<"\n";
  
  IntegerVector both_child_term(internal_nodes_prop.size());
 // std::cout<<"Internal_nodes_prop"<<internal_nodes_prop.size();
  /*for(int  i = 0; i < internal_nodes_prop.size(); ++i){
   std::cout <<"in for loop"<< internal_nodes_prop[i] << ' ';
  }*/
  for(int k=0;k<internal_nodes_prop.size();k++){
   // std::cout<<"Row index is"<<" "<<k<<" "<<tree_table(internal_nodes_prop[k]-1,1);
    //std::cout<<"tree index is"<<" "<<k<<" "<<tree_table(tree_table(internal_nodes_prop[k]-1,1)-1,4)<<"\n";
   // std::cout<<"tree index is"<<tree_table(tree_table(internal_nodes_prop[k]-1,1),4)<<"\n";
  // std::cout<<"logical test"<<(tree_table(tree_table(internal_nodes_prop[k]-1,0),4)==-1)<<"\n";
    if((tree_table(tree_table(internal_nodes_prop[k]-1,0)-1,4)==
    tree_table(tree_table(internal_nodes_prop[k]-1,1)-1,4)) && (tree_table(tree_table(internal_nodes_prop[k]-1,0),4)==-1)){
      //std::cout<<"Both children are terminal! \n";
      both_child_term[k]=1;
    
    }else{
      both_child_term[k]=0;
    }
  }
 // std::cout<<"both_child_term"<<both_child_term[3]<<"\n";
  IntegerVector nodes_to_prune;
  
  for(int k=0;k<both_child_term.size();k++){
    if(both_child_term[k]==1){
      nodes_to_prune.push_back(internal_nodes_prop[k]);     
    }
  }  
 // std::cout<<"both_child_term"<<both_child_term[3]<<"\n";
  return(nodes_to_prune);
}