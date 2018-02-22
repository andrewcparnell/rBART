//########################################################################################//
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
NumericVector get_grow_obs(NumericMatrix xmat,NumericVector grow_obs,int split_var){
  NumericVector sv_col=xmat(_,split_var-1);
  NumericVector get_min;
  
  for(int i=0;i<grow_obs.size();i++){    
    get_min.push_back(sv_col[grow_obs[i]]);
  }
  
return(get_min);
}



#include <Rcpp.h>
using namespace Rcpp ;
// [[Rcpp::export]]
NumericMatrix get_subset(NumericMatrix xmat,NumericVector grow_obs){
  NumericMatrix B(grow_obs.size(),xmat.ncol());
  for(int i=0;i<grow_obs.size();i++){
    B(i,_)=xmat(i,_);
  }
  
return(B);
}

//******************************************************************************************/
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
    
    if(sv_col[obs_to_update[i]]<=split_point){
     
      ld_obs.push_back(obs_to_update[i]);
  
    }
    if(sv_col[obs_to_update[i]]>=split_point){  
     
      rd_obs.push_back(obs_to_update[i]);

    }
  }
  
  daughter_obs[0]=ld_obs;
  daughter_obs[1]=rd_obs;
  return(daughter_obs);
}

//#######################################################################################//
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix add_rows(NumericMatrix prior_tree_table_temp,int grow_node){
  arma::mat M=Rcpp::as<arma::mat>(prior_tree_table_temp);
  M(grow_node-1,5)=0;
  M(grow_node-1,6)=0;
  M(grow_node-1,0)=grow_node+1;
  M(grow_node-1,1)=grow_node+2;
  M.insert_rows(grow_node,2);
  M(grow_node,4)=-1;
  M(grow_node+1,4)=-1;
  NumericMatrix t=as<NumericMatrix>(wrap(M));
  IntegerVector rname=seq_len(t.nrow());
  
   List dimnms = // two vec. with static names
    List::create(rname,
                 CharacterVector::create("left daughter","right daughter","split var","split point","status","mean","std dev"));
    // and assign it
    t.attr("dimnames") = dimnms;
  return(t);
}

//#######################################################################################//
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
IntegerMatrix addcol(IntegerMatrix prior_tree_matrix_temp,int grow_node,NumericVector ld_obs,NumericVector rd_obs){
  int ncol=prior_tree_matrix_temp.ncol();
  arma::mat M=Rcpp::as<arma::mat>(prior_tree_matrix_temp);
  M.insert_cols(ncol,1);
  for(int i =0;i<ld_obs.size();i++){
    M(ld_obs[i],ncol)=grow_node+1;
  }
  for(int i =0;i<rd_obs.size();i++){
    M(rd_obs[i],ncol)=grow_node+2;
  }

  IntegerMatrix t=as<IntegerMatrix>(wrap(M));
  return(t);
} 


#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

NumericMatrix set_daughter_to_end_tree(int grow_node,NumericMatrix prior_tree_table_temp,double left_daughter){
  int nrow=prior_tree_table_temp.nrow();
  arma::mat M=Rcpp::as<arma::mat>(prior_tree_table_temp);
  M(grow_node-1,5)=0;
  M(grow_node-1,6)=0;
  M.insert_rows(nrow,2);

  M(grow_node-1,0)=left_daughter;
  M(grow_node-1,1)=left_daughter+1;
  M(left_daughter-1,4)=-1;
  M(left_daughter,4)=-1;

  NumericMatrix s=as<NumericMatrix>(wrap(M));
  IntegerVector rname=seq_len(s.nrow());
  
   List dimnms = // two vec. with static names
    List::create(rname,
                 CharacterVector::create("left daughter","right daughter","split var","split point","status","mean","std dev"));
    // and assign it
    s.attr("dimnames") = dimnms;

  return(s);
}
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

IntegerMatrix set_daughter_to_end_mat(double d,IntegerMatrix prior_tree_matrix_temp,double left_daughter,NumericVector ld_obs,NumericVector rd_obs){
  int ncol_mat=prior_tree_matrix_temp.ncol();
  arma::mat N=Rcpp::as<arma::mat>(prior_tree_matrix_temp);
         
  if(d+1==ncol_mat){
  //  #now need to update prior_tree_matrix
 //   #first need to insert extra column for the new split node

    N.insert_cols(ncol_mat,1); 
 
    for(int i =0;i<ld_obs.size();i++){
      N(ld_obs[i],d+1)=left_daughter;
     
    }
    for(int i =0;i<rd_obs.size();i++){
      N(rd_obs[i],d+1)=left_daughter +1;
    }
  }else{
      for(int i =0;i<ld_obs.size();i++){
      N(ld_obs[i],d+1)=left_daughter;
    }
    for(int i =0;i<rd_obs.size();i++){
      N(rd_obs[i],d+1)=left_daughter +1;
    }  
  }

  IntegerMatrix t=as<IntegerMatrix>(wrap(N));
  
  return(t);
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
IntegerVector remove_zero(IntegerVector nodes_at_depth){
  IntegerVector ret;
  for(int i=0;i<nodes_at_depth.size();i++){
    if(nodes_at_depth[i]!=0){
      ret.push_back(nodes_at_depth[i]);
    }
  }
  return(ret);
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
IntegerVector get_gnp(IntegerVector nodes_at_depth,int grow_node){
 
  IntegerVector grow_node_pos;
  for(int i =0;i<nodes_at_depth.size();i++){
    if(nodes_at_depth[i]==grow_node){
      grow_node_pos.push_back(i);
    } 
  }
  return(grow_node_pos);  
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

int find_prev_nonterm(IntegerVector find_nonterm,IntegerVector prev){
  
  int ret=0;
  for(int j=0;j<prev.size();j++){
    for(int i=0;i<find_nonterm.size();i++){
      if(prev[j]==find_nonterm[i]){
        ret+=1;
      }
    }
  }  
  return(ret);
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

NumericVector find_nodes_to_update(NumericVector all_ld,double left_daughter){
  NumericVector gr_ld;
  for(int i=0;i<all_ld.size();i++){
    if(all_ld[i]>=left_daughter){
      gr_ld.push_back(i);
    }
  }
  return(gr_ld);
}

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix set_tree_to_middle(NumericVector node_to_update,NumericMatrix prior_tree_table_temp,int grow_node,double left_daughter){
  for(int i=0;i<node_to_update.size();i++){
    if(prior_tree_table_temp(node_to_update[i],0) && prior_tree_table_temp(node_to_update[i],1)!=0){
      prior_tree_table_temp(node_to_update[i],0)+=2;
      prior_tree_table_temp(node_to_update[i],1)+=2;
    }
  }
 
      prior_tree_table_temp(grow_node-1,5)=0;
      prior_tree_table_temp(grow_node-1,6)=0;
      
      arma::mat M=Rcpp::as<arma::mat>(prior_tree_table_temp);
      M.insert_rows(left_daughter-1,2);
      M(left_daughter-1,4)=-1;
      M(left_daughter,4)=-1;      
  
      M(grow_node-1,0)=left_daughter;
      M(grow_node-1,1)=left_daughter+1;
      NumericMatrix t=as<NumericMatrix>(wrap(M));
      IntegerVector rname=seq_len(t.nrow());
  
      List dimnms = // two vec. with static names
      List::create(rname,
                 CharacterVector::create("left daughter","right daughter","split var","split point","status","mean","std dev"));
    // and assign it
    t.attr("dimnames") = dimnms;
  return(t);
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
IntegerMatrix find_obs_to_update_grow(IntegerMatrix prior_tree_matrix_temp,double left_daughter,double d,NumericVector ld_obs,NumericVector rd_obs){
  IntegerVector rows_obs;
  IntegerVector cols_obs;
  
  for(int i=0;i<prior_tree_matrix_temp.nrow();i++){
    for(int j=0;j<prior_tree_matrix_temp.ncol();j++){
      if(prior_tree_matrix_temp(i,j)>=left_daughter){
        rows_obs.push_back(i);
        cols_obs.push_back(j);
      }
    }
  }
  
  if(rows_obs.size()!=0){
    for(int k=0;k< rows_obs.size();k++){
      if(prior_tree_matrix_temp(rows_obs[k],cols_obs[k])<left_daughter){
            prior_tree_matrix_temp(rows_obs[k],cols_obs[k])=prior_tree_matrix_temp(rows_obs[k],cols_obs[k]);
      }else if(prior_tree_matrix_temp(rows_obs[k],cols_obs[k])==0){
          prior_tree_matrix_temp(rows_obs[k],cols_obs[k])=0;
       }else{   
         int temp=prior_tree_matrix_temp(rows_obs[k],cols_obs[k])+2;
          prior_tree_matrix_temp(rows_obs[k],cols_obs[k])=temp;
        }
    }
  }
   // #now need to update prior_tree_matrix
  // #first need to insert extra column for the new split node
  if(prior_tree_matrix_temp.ncol()>d+1){
    arma::mat M=Rcpp::as<arma::mat>(prior_tree_matrix_temp);
    M.insert_cols(prior_tree_matrix_temp.ncol(),1);  
    IntegerMatrix prior_tree_matrix=as<IntegerMatrix>(wrap(M));
  }
  for(int i =0;i<ld_obs.size();i++){
    
      prior_tree_matrix_temp(ld_obs[i],d+1)=left_daughter;
     
    }
    for(int i =0;i<rd_obs.size();i++){
       prior_tree_matrix_temp(rd_obs[i],d+1)=left_daughter+1;  
    }    
 
  return(prior_tree_matrix_temp);
}
//#######################################################################################//
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

List grow_tree(NumericMatrix xmat,NumericVector y,IntegerMatrix prior_tree_matrix,int grow_node,NumericMatrix prior_tree_table)
{
  IntegerMatrix prior_tree_matrix_temp=clone(prior_tree_matrix);
  NumericMatrix prior_tree_table_temp=clone(prior_tree_table);
  double yy=xmat.ncol();
  IntegerVector xx=seq_len(yy);
  
  NumericVector terminal_nodes=find_term_nodes(prior_tree_table_temp);
  IntegerVector splitvar1=csample_num(xx,1,false);
  int splitvar=splitvar1[0];
  
  // now get the min and max values of that variable for the observations within the chosen 
  //  terminal node 
  
 
  NumericVector grow_obs=find_term_obs(prior_tree_matrix_temp,grow_node);
  
  NumericVector get_min=get_grow_obs(xmat,grow_obs,splitvar);
  //depth of tree at current terminal node
  NumericVector d1=unique(find_term_cols(prior_tree_matrix_temp,grow_node));
  double d=d1[0];
  double min_point=min(get_min);
  double max_point=max(get_min);
  //choose a random split point in the range of the split variable for the observations in 
  //the chosen terminal node
  double splitpoint=R::runif(min_point,max_point);
//List add_daughter(2);
  //#get the subset of data available for the terminal node we are growing:
  
  NumericMatrix data_curr_node=get_subset(xmat,grow_obs);
  
  prior_tree_table_temp(grow_node-1,3)=splitpoint;
  prior_tree_table_temp(grow_node-1,2)=splitvar;
  prior_tree_table_temp(grow_node-1,4)=1;
  
  //  #need the mean and standard deviation for y are already reported in prior_tree_table
 //  #get data subset for left and right daughter nodes
  List daughter_obs=get_daughter_obs(xmat,grow_obs,splitvar,splitpoint);
  NumericVector ld_obs=daughter_obs[0];
  NumericVector rd_obs=daughter_obs[1];
  
  if(prior_tree_table_temp.nrow()==grow_node){
    
  //if the grow node is the final node in the tree:
    
  //if the terminal node chosen is the last one in the tree
    prior_tree_table_temp=add_rows(prior_tree_table_temp,grow_node);
   
   // #now need to add new column to the prior_tree_matrix matrix
    prior_tree_matrix_temp=addcol(prior_tree_matrix_temp,grow_node,ld_obs,rd_obs);  
    
  }else{
    
  //  #if grow node is in the middle of the tree
    IntegerVector nodes_d;
    nodes_d=prior_tree_matrix_temp(_,d);
    IntegerVector nodes_at_depth=sort_unique(nodes_d);
    //nodes_at_depth<-sort(nodes_at_depth)
    IntegerVector nodes_at_depth1=remove_zero(nodes_at_depth);
    IntegerVector gn_pos=get_gnp(nodes_at_depth1, grow_node);
   
    IntegerVector prev;
    for(int i=0;i<gn_pos.size();i++){
       prev.push_back(nodes_at_depth1[i]);
    }
    IntegerVector find_nonterm=find_internal_nodes(prior_tree_table);
    
    int prev_nonterm=find_prev_nonterm(find_nonterm,prev);
   
    // #node number for left daughter will be grow_node + 2 for each previous node that was split and +1 for each
    //#node after current grow node at depth d
    double left_daughter=grow_node +2*(prev_nonterm)+(nodes_at_depth1.size()-gn_pos[0]);
    
    NumericVector node_to_update=find_nodes_to_update(prior_tree_table(_,1),left_daughter);
   
    //  #in the tree matrix need to increase the node number of nodes after the grow node by two (because we added in 2 daughter nodes to grow node)
    //   #do this for all observations except those that already belong to a terminal node (a value of 0)
    //   #  prior_tree_matrix[obs_to_update[,1],obs_to_update[,2]]<-ifelse(prior_tree_matrix[obs_to_update[,1],obs_to_update[,2]]==0,0,prior_tree_matrix[obs_to_update[,1],obs_to_update[,2]]+2))
    if(node_to_update.size()==0){
      
      //#daughter nodes need to be added to the end of the table not in the center of it
      //add_daughter=set_daughter_to_end(grow_node,d,prior_tree_table_temp,prior_tree_matrix_temp,left_daughter,ld_obs,rd_obs);
     prior_tree_table_temp=set_daughter_to_end_tree(grow_node,prior_tree_table_temp,left_daughter);
     prior_tree_matrix_temp=set_daughter_to_end_mat(d,prior_tree_matrix_temp,left_daughter,ld_obs,rd_obs);
    }else{
   //   #if the daughter node number already exists in the tree and existing node numbers have to
   //   #be updated
      prior_tree_table_temp=set_tree_to_middle(node_to_update,prior_tree_table_temp,grow_node,left_daughter);
      prior_tree_matrix_temp=find_obs_to_update_grow(prior_tree_matrix_temp,left_daughter,d,ld_obs,rd_obs);    
      
    }
  }
  List ret(2);
  ret[0]=prior_tree_table_temp;
  ret[1]=prior_tree_matrix_temp;
  
  return(ret);
}