#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

List get_daughter_obs(NumericMatrix xmat,NumericVector obs_to_update,int split_var,double split_point){
  IntegerVector ld_obs;
  IntegerVector rd_obs;
  NumericVector sv_col;
  List daughter_obs(2);
  //NumericVector sv_obs;
  //sv_obs=xmat(obs_to_update,_);
  sv_col=xmat(_,split_var-1);
 
  
  for(int i=0;i<obs_to_update.size();i++){
    std::cout<<sv_col[obs_to_update[i]-1]<<"\n";
    
    if(sv_col[obs_to_update[i]-1]<=split_point){
      std::cout<<"xmat col2"<<sv_col[obs_to_update[i]-1]<<"\n";
     
      ld_obs.push_back(obs_to_update[i]-1);
      //std::cout<<"ld_obs"<<ld_obs[0]<<"\n";
    }
    if(sv_col[obs_to_update[i]-1]>=split_point){
      std::cout<<"xmat col2"<<sv_col[obs_to_update[i]-1]<<"\n";
     
      rd_obs.push_back(obs_to_update[i]-1);
      //std::cout<<"ld_obs"<<ld_obs[0]<<"\n";
    }
  }
  //ld_obs=xmat(obs_to_update[,1],)
  //ld_obs<-ld_obs[which(ld_obs[,split_var]<=split_point),]
  daughter_obs[0]=ld_obs;
  daughter_obs[1]=rd_obs;
  return(daughter_obs);
}

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


#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

NumericVector find_term_cols(IntegerMatrix tree_matrix,int terminal_node){
  NumericVector term_cols;
for(int i=0;i<tree_matrix.nrow();i++){
      for(int j=0;j<tree_matrix.ncol();j++){
        
        if(tree_matrix(i,j)==terminal_node){
          term_cols.push_back(j);  
        }
      }
    }
    return(term_cols);
}

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
//***************************************************************************************/
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

IntegerMatrix update_tree_matrix(NumericMatrix xmat,int swap_parent,NumericMatrix treedata,IntegerMatrix tree_matrix){
  int left_daughter;
  int right_daughter;
  left_daughter=treedata(swap_parent-1,0);
  right_daughter=treedata(swap_parent-1,1);
  int ind_nodes_to_update=1;
  IntegerVector nodes_to_update;
  nodes_to_update[0]=swap_parent;  
  int split_var;
  double split_point;
  //while(ind_nodes_to_update!=0){
  //go in here when there is still a node that was split so observations in children 
  //nodes need to be reassigned
    for(int i=0;i<ind_nodes_to_update;i++){
    
      NumericVector obs_to_update=find_term_obs(tree_matrix,nodes_to_update[i]);
    
      NumericVector term_cols=unique(find_term_cols(tree_matrix,nodes_to_update[i]));
     // tree_matrix(obs_to_update,unique(find_term_cols(tree_matrix,nodes_to_update[i]))+1)=999;

      std::cout<<"obs_to_update"<<obs_to_update.size()<<"cols to update"<<term_cols[0]<<"\n";
      std::cout<<"left_daughter"<<left_daughter<<"right_daughter"<<right_daughter<<"\n";
      
      for(int j =0;j<obs_to_update.size();j++){
        tree_matrix(obs_to_update[j],term_cols[0]+1)=999;
      }
      left_daughter=treedata(nodes_to_update[i]-1,0);    
      right_daughter=left_daughter+1;
      split_var=treedata(nodes_to_update[i]-1,2);
      split_point=treedata(nodes_to_update[i]-1,3);
      std::cout<<"ld "<<left_daughter<<"rd "<<right_daughter<<"sv "<<split_var<<"sp "<<split_point<<"\n";  
       we haven't got to the end of the branch yet
        //follow it down and keep refining the subset until we get to a terminal node
      IntegerVector ld_obs;
      IntegerVector rd_obs;
      
      ld_obs=get_daughter_obs(xmat,obs_to_update,split_var,split_point)[0];
      rd_obs=get_daughter_obs(xmat,obs_to_update,split_var,split_point)[1];
      std::cout<<"ld "<<ld_obs.size()<<"rd "<<rd_obs.size()<<"\n"; 
     
    }
  return(tree_matrix);
}