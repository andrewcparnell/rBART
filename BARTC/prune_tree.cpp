#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

// [[Rcpp::export]]
IntegerVector csample_num( IntegerVector x,
                           int size,
                           bool replace,
                           NumericVector prob = NumericVector::create()
                           ) {
  RNGScope scope;
  IntegerVector ret = RcppArmadillo::sample(x, size, replace, prob);
  return ret;
}

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List shed(NumericMatrix prior_tree_table_temp,IntegerMatrix prior_tree_matrix_temp){
  arma::mat M=Rcpp::as<arma::mat>(prior_tree_table_temp);
  M.shed_rows(1,2); 
  M(0,4)=-1;
  M(0,0)=0;
  M(0,1)=0;
  M(0,2)=0;
  M(0,3)=0;
  arma::mat N=Rcpp::as<arma::mat>(prior_tree_matrix_temp);
  N.shed_col(1);
  NumericMatrix s=as<NumericMatrix>(wrap(N));
  NumericMatrix t=as<NumericMatrix>(wrap(M));
  IntegerVector rname=seq_len(t.nrow());
  
  List dimnms = // two vec. with static names
  List::create(rname,
         CharacterVector::create("left daughter","right daughter","split var","split point","status","mean","std dev"));
  // and assign it
  t.attr("dimnames") = dimnms;
  
  List ret(2);
  ret[0]=t;
  ret[1]=s;
  return(ret);
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

IntegerVector find_internal_nodes(NumericMatrix treetable){
  IntegerVector internal_nodes;

  for(int l=0;l<treetable.nrow();l++){
    
    if(treetable(l,4)==1){
      internal_nodes.push_back(l+1);
    }
  }
  IntegerVector internal_nodes_sort = clone(internal_nodes);
  std::sort(internal_nodes.begin(), internal_nodes.end());
  
  
  return(internal_nodes_sort);
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

IntegerVector find_prune_node(IntegerVector internal_nodes,NumericMatrix prior_tree_table_temp){
    IntegerVector nodes_to_prune;
        
    for(int i=0;i<internal_nodes.size();i++)
    {
      //look at the tree table to see if the left and right daughter of each internal node
      //are terminal. If they are then the internal node is prunable.
      
      if((prior_tree_table_temp(prior_tree_table_temp(internal_nodes[i]-1,0)-1,4))==(prior_tree_table_temp(prior_tree_table_temp(internal_nodes[i]-1,1)-1,4))&& (prior_tree_table_temp(prior_tree_table_temp(internal_nodes[i]-1,0)-1,4)==-1)){
        nodes_to_prune.push_back(internal_nodes[i]);
      }
    }
  return(nodes_to_prune);
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

NumericVector find_term_cols(IntegerMatrix tree_matrix_temp,int terminal_node){
  NumericVector term_cols;
for(int i=0;i<tree_matrix_temp.nrow();i++){
      for(int j=0;j<tree_matrix_temp.ncol();j++){
        
        if(tree_matrix_temp(i,j)==terminal_node){
          term_cols.push_back(j);  
        }
      }
    }
    return(term_cols);
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

NumericMatrix set_prune_table(NumericMatrix prior_tree_table_temp,int prune_node,int ld,int rd){
  
  prior_tree_table_temp(prune_node-1,0)=0;
  prior_tree_table_temp(prune_node-1,1)=0;
  prior_tree_table_temp(prune_node-1,2)=0;
  prior_tree_table_temp(prune_node-1,3)=0;
  prior_tree_table_temp(prune_node-1,4)=-1;
  return(prior_tree_table_temp);
}   

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

NumericMatrix set_prune_tree_end(NumericMatrix prior_tree_table_temp,int ld,int rd){
  arma::mat M=Rcpp::as<arma::mat>(prior_tree_table_temp);

  M.shed_rows(ld-1,rd-1);
  NumericMatrix s=as<NumericMatrix>(wrap(M));  
  IntegerVector rname=seq_len(s.nrow());
   List dimnms = // two vec. with static names
    List::create(rname,
                 CharacterVector::create("left daughter","right daughter","split var","split point","status","mean","std dev"));
    // and assign it
    s.attr("dimnames") = dimnms;
  return(s);
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

NumericVector find_term_obs(IntegerMatrix tree_matrix_temp,int terminal_node){
  NumericVector term_obs;
  for(int i=0;i<tree_matrix_temp.nrow();i++){
      for(int j=0;j<tree_matrix_temp.ncol();j++){
        
        if(tree_matrix_temp(i,j)==terminal_node){
          term_obs.push_back(i);  
        }
      }
    }
    return(term_obs);
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

IntegerMatrix set_prune_mat_end(IntegerMatrix prior_tree_matrix_temp,int ld,int rd,double d){
  NumericVector ld_obs=find_term_obs(prior_tree_matrix_temp,ld);
  NumericVector ld_cols=find_term_cols(prior_tree_matrix_temp,ld);
  NumericVector rd_obs=find_term_obs(prior_tree_matrix_temp,rd);
  NumericVector rd_cols=find_term_cols(prior_tree_matrix_temp,rd);
  
  for(int i=0;i<ld_obs.size();i++){
    prior_tree_matrix_temp(ld_obs[i],ld_cols[i])=0;
  }
  for(int i=0;i<rd_obs.size();i++){
    prior_tree_matrix_temp(rd_obs[i],rd_cols[i])=0;
  }
  
  if(sum(prior_tree_matrix_temp(_,d+1))==0 ){
  
  // #if the nodes to be pruned were the only ones at that depth then we should delete this depth from tree
    arma::mat M=Rcpp::as<arma::mat>(prior_tree_matrix_temp);
    M.shed_col(d+1);
    prior_tree_matrix_temp=as<IntegerMatrix>(wrap(M));
    
  }
  
  return(prior_tree_matrix_temp);
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix set_prune_table_middle(NumericMatrix prior_tree_table_temp,int ld,int rd){
  IntegerVector node_to_update;
  for(int i=0;i<prior_tree_table_temp.nrow();i++){
    if(prior_tree_table_temp(i,0)>ld){
      node_to_update.push_back(i);
    }
  }
  //for each node to update, change their node assignment
  for(int i=0;i<node_to_update.size();i++){
    
    prior_tree_table_temp(node_to_update[i],0)=prior_tree_table_temp(node_to_update[i],0)-2;
    prior_tree_table_temp(node_to_update[i],1)=prior_tree_table_temp(node_to_update[i],1)-2;
    
    arma::mat M=Rcpp::as<arma::mat>(prior_tree_table_temp);
    M.shed_rows(ld-1,rd-1);  
    prior_tree_table_temp=as<NumericMatrix>(wrap(M));
    IntegerVector rname=seq_len(prior_tree_table_temp.nrow());
   
    List dimnms = // two vec. with static names
    List::create(rname,
                 CharacterVector::create("left daughter","right daughter","split var","split point","status","mean","std dev"));
    // and assign it
    prior_tree_table_temp.attr("dimnames") = dimnms;
  }
 
  return(prior_tree_table_temp);     
}


#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
IntegerMatrix set_prune_mat_middle(IntegerMatrix prior_tree_matrix_temp,int ld,int rd){
  IntegerVector obs_to_update;
  IntegerVector cols_to_update;
  NumericVector ld_obs=find_term_obs(prior_tree_matrix_temp,ld);
  NumericVector ld_cols=find_term_cols(prior_tree_matrix_temp,ld);
  NumericVector rd_obs=find_term_obs(prior_tree_matrix_temp,rd);
  NumericVector rd_cols=find_term_cols(prior_tree_matrix_temp,rd);
  
  for(int i=0;i<prior_tree_matrix_temp.nrow();i++){
    for(int j=0;j<prior_tree_matrix_temp.ncol();j++){
      if(prior_tree_matrix_temp(i,j)>rd){
        obs_to_update.push_back(i);
        cols_to_update.push_back(j);
      }
    }  
  }   

  //#in the tree matrix need to increase the node number of nodes after the grow node by two (because we added in 2 daughter nodes to grow node)
  // #do this for all observations except those that already belong to a terminal node (a value of 0)
  for(int i=0;i<ld_obs.size();i++){
    prior_tree_matrix_temp(ld_obs[i],ld_cols[i])=0;
  }
  for(int i=0;i<rd_obs.size();i++){
    prior_tree_matrix_temp(rd_obs[i],rd_cols[i])=0;
  }
      
  if(obs_to_update.size()!=0){
    for(int k=0; k<obs_to_update.size();k++){
      if((prior_tree_matrix_temp(obs_to_update[k],cols_to_update[k]))==0){
        prior_tree_matrix_temp(obs_to_update[k],cols_to_update[k])=0;
      }else{
        prior_tree_matrix_temp(obs_to_update[k],cols_to_update[k])=prior_tree_matrix_temp(obs_to_update[k],cols_to_update[k])-2;
      }
    }
  }
  return(prior_tree_matrix_temp); 
}
/*########################################################################################*/
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

List prune_tree(NumericMatrix xmat,IntegerMatrix prior_tree_matrix,NumericMatrix prior_tree_table){
  IntegerMatrix prior_tree_matrix_temp=clone(prior_tree_matrix);
  NumericMatrix prior_tree_table_temp=clone(prior_tree_table);
  
  List ret(2);

  if(prior_tree_table_temp.nrow()==1){
   
   ret[0]=prior_tree_table_temp;
   ret[1]=prior_tree_matrix_temp;
  
  }else if(prior_tree_table_temp.nrow()==3){
      
       ret=shed(prior_tree_table_temp, prior_tree_matrix_temp);
   
    }else{
      IntegerVector internal_nodes=find_internal_nodes(prior_tree_table_temp);
      
      //both_child_term<-rep(0,length(internal_nodes))
      IntegerVector nodes_to_prune=find_prune_node(internal_nodes,prior_tree_table_temp);
      
      
      IntegerVector xx=seq_len(nodes_to_prune.size());
     
      IntegerVector  choose_node=csample_num(xx,1,false);
      int choose_node1=choose_node[0]; 
      int prune_node=nodes_to_prune[choose_node1-1];
      
      int ld=prior_tree_table_temp(prune_node-1,0);
      int rd=prior_tree_table_temp(prune_node-1,1);
      
    // #get depth of tree at current prune node
      NumericVector d1=unique(find_term_cols(prior_tree_matrix_temp,prune_node));
      double d=d1[0];
     
    //#get the subset of data available for the terminal node we are growing:
    //data_curr_node<-xmat[which(prior_tree_matrix[,d]==prune_node),]
      prior_tree_table_temp=set_prune_table(prior_tree_table_temp,prune_node,ld,rd);
       
      if(prior_tree_table_temp.nrow()==rd){
        
      //#if the prune node is the final internal node in the tree:
        prior_tree_table_temp=set_prune_tree_end(prior_tree_table_temp,ld,rd);
        prior_tree_matrix_temp=set_prune_mat_end(prior_tree_matrix_temp,ld,rd,d);
      }else{
      //#if the two pruned nodes are not the last two nodes in the tree have to update the node numbers of subsequent nodes
      prior_tree_table_temp=set_prune_table_middle(prior_tree_table_temp,ld,rd);
      prior_tree_matrix_temp=set_prune_mat_middle(prior_tree_matrix_temp,ld,rd);
      
      }  
      //#internal nodes where LD is terminal  
      ret[0]=prior_tree_table_temp;
      ret[1]=prior_tree_matrix_temp;
    }
  return(ret);
}  