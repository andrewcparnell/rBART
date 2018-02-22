#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

double get_tree_prior(NumericMatrix tree_table,IntegerMatrix tree_matrix,double alpha,double beta){
  double propsplit=1;
  IntegerVector d;
  IntegerVector int_nodes_index;
  NumericVector internal_nodes_prop;
  
  for(int i=0;i<tree_table.nrow();i++){
   // std::cout<<"i is"<<" "<<i<<"node status is"<<node_status[i]<<"\n";
    if(tree_table(i,4)==1){
    //std::cout<<"NODE STATUS=1!!"<<"\n";
      internal_nodes_prop.push_back(i+1); 
    }
  }
  
  for(int k=0;k<internal_nodes_prop.size();k++){ 
   // std::cout<<"internal nores are"<<internal_nodes_prop[k]<<"\n";
    for(int i=0;i<tree_matrix.nrow();i++){
      for(int j=0;j<tree_matrix.ncol();j++){
        
        if(tree_matrix(i,j)==internal_nodes_prop[k]){
          int_nodes_index.push_back(j+1);  
        }
      }
    }
   // std::cout<<"internal node index is"<<int_nodes_index.size()<<"\n";
    
    if(int_nodes_index.size()!=0){
      d=unique(int_nodes_index);
      double d1=d[0];
   //   std::cout<<"Node depth is"<<d[0]<<"\n";
      //double prior_p=alpha*(1+d1);
     // std::cout<<"prior p is"<<prior_p <<"beta is"<<beta<<"\n";
      propsplit*=alpha*pow((d1+1),-beta) ;
      //propsplit=1.0/propsplit;
    //  std::cout<<"propsplit "<<k<<" time around is"<<propsplit<<"\n";
    }
    IntegerVector temp;
    int_nodes_index=temp;
  } 
  return(propsplit);
}