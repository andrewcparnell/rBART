#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix remove_curr_col(NumericMatrix predy,int i){
  arma::mat M=Rcpp::as<arma::mat>(predy);
  M.shed_col(i);
  NumericMatrix s=as<NumericMatrix>(wrap(M));
return(s);
}
//####################################################################################//
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
//##########################################################################################//
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
//############################################################################################//
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

//############################################################################################//
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
//############################################################################################//
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
//############################################################################################//
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
//##########################################################################################//
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
//##########################################################################################//
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

NumericVector find_term_nodes(NumericMatrix tree_table){
  NumericVector terminal_nodes;
  
  for(int i=0;i<tree_table.nrow();i++){
    if(tree_table(i,4)==-1){
      terminal_nodes.push_back(i+1);
    }
  }
  return(terminal_nodes);
}  

//#######################################################################################//
//##########################################################################################//
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
//########################################################################################//

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
//############################################################################################//
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
//############################################################################################//
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
//############################################################################################//
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
//###########################################################################################//
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

List get_daughter_obs(NumericMatrix xmat,NumericVector obs_to_update,int split_var,double split_point){
  IntegerVector ld_obs;
  IntegerVector rd_obs;
  NumericVector sv_col;
  List daughter_obs(2);
  //NumericVector sv_obs;
 
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
//############################################################################################//
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
//##########################################################################################//
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
//##########################################################################################//
//##########################################################################################//
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
//#########################################################################################//
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
//############################################################################################//

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
//##########################################################################################//


//############################################################################################//
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

//******************************************************************************************//
//get the log likelihood
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

double likelihood_function(NumericVector y_temp,NumericMatrix treetable_temp,IntegerMatrix obs_to_nodes_temp,double a,double mu,double nu,double lambda){
  
  IntegerVector terminal_nodes;

  for(int l=0;l<treetable_temp.nrow();l++){
    
    if(treetable_temp(l,4)==-1){
      terminal_nodes.push_back(l+1);
    }
  }
  double b=terminal_nodes.size();
  //NumericVector n(b);
  NumericVector sumr(b);
  NumericVector sumsqr(b);
  NumericVector term_obs;
  int ni;
  IntegerVector n(b);
  NumericVector t(b);
  NumericVector s(b);
 
  //parameters that do change for each node
  //terminal_nodes
  double p2=0;
  double p3=0;
  double p31=0;
  double p32=0;

  for(int i=0;i< b;i++){
    //number of observations in terminal nodes
    term_obs=find_term_obs(obs_to_nodes_temp,terminal_nodes[i]);
    ni=term_obs.size();
    n[i]=ni;
    
    double temp1=0;
    double temp2=0;
   
    
    for(int j=0;j<ni;j++){
      temp1+=y_temp[term_obs[j]];    
      temp2+=pow(y_temp[term_obs[j]],2);    
    }

   sumr[i]=temp1;
    sumsqr[i]=temp2;
    if(ni!=0){
    s[i] = (sumsqr[i]-pow(sumr[i],2)/n[i]);
     //s=s[!is.na(s)]
    t[i] = (n[i]*a/(n[i]+a))*pow((sumr[i]/n[i]-mu),2);
     //t=t[!is.na(t)]
    }else{
      s[i] =0;
      t[i] =0;
    }
 
  p2+=0.5*log(n[i]+a);
  p31+=(n[i]);
  p32+=s[i]+t[i]+nu*lambda;
    
  } 
  
  double la=log(a);
  double p1=(b/2)*la;
  p3=(p31+nu)/2*log(p32);
 
  double tree_log_lik=p1 - p2 -p3;
  //double tree_log_lik=(b/2)*log(a) - sum(0.5*log(n+a)) -(sum(n)+nu)/2*log(sum((s+t)+nu*lambda));

 //if the chosen change in split points gives a tree with an empty node then
  //multiply the likelihood by a big number so this tree won't be chosen
  for(int i=0;i<n.size();i++){
    if(n[i]<=5){
      tree_log_lik=tree_log_lik-(100000);
    }
    else{
      tree_log_lik=tree_log_lik;
    }
  }
return(tree_log_lik);
}
//############################################################################################//
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix start_predy(int num_trees,int n,double start_mean,double start_sd){
  NumericMatrix predy(n,num_trees);
  int ncol=predy.ncol();
  int nrow=predy.nrow();
  for(int i=0;i<ncol;i++){
    for(int j=0;j<nrow;j++){
      predy(j,i)=R::rnorm(start_mean,start_sd);
    }  
  }
  return(predy);
}
//#############################################################################################//
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector scale_response(double a,double b,double c,double d,NumericVector y){
  NumericVector y_scaled = -((-b*c+a*d)/(-a+b))+((-c+d)*y/(-a+b));
  return(y_scaled);
}
//#############################################################################################//
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List start_tree(int num_trees,double start_mean,double start_sd){
  List pri_tree(num_trees);
  NumericMatrix treemat(1,7);
  
  
  for(int i =0;i<num_trees;i++){
    double rand=R::rnorm(start_mean,start_sd);
    NumericVector testrow = NumericVector::create(0,0,0,0,-1,rand,0);
    for(int k=0;k<1;k++){
      for(int j=0;j<7;j++){
        treemat(k,j)=testrow[j];
      }
  }
    //testrow[5]=rand;
    List dimnms = // two vec. with static names
    List::create(CharacterVector::create("1"),
                 CharacterVector::create("left daughter","right daughter","split var","split point","status","mean","std dev"));
    // and assign it
    treemat.attr("dimnames") = dimnms;
    pri_tree[i]=treemat;
  }
                   
  return(pri_tree);
}
//#############################################################################################//
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
//#############################################################################################//
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
//############################################################################################//
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
//############################################################################################//
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
//##########################################################################################//
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

//##########################################################################################//
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
//##########################################################################################//
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
IntegerMatrix set_daughter(int left_daughter,int right_daughter,IntegerVector ld_obs,IntegerVector rd_obs,IntegerMatrix tree_matrix_temp,double term_cols){
  for(int i=0;i<ld_obs.size();i++){
    tree_matrix_temp(ld_obs[i],term_cols+1)=left_daughter;    
  }
  for(int i=0;i<rd_obs.size();i++){
    tree_matrix_temp(rd_obs[i],term_cols+1)=right_daughter;    
  }    
  return(tree_matrix_temp);
}
//############################################################################################//
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List update_changetree_matrix(NumericMatrix xmat,NumericMatrix treedata_temp,IntegerMatrix tree_matrix_temp,double split_var,double split_point,int change_node){
  List ret(2);
  int left_daughter;
  int right_daughter;
  IntegerVector ld_obs;
  IntegerVector rd_obs;
  int ind_nodes_to_update=1;
  int temp_ind_nodes_to_update=ind_nodes_to_update;
  IntegerVector temp_nodes_to_update;
 //#assign new split var and point to the treedata
  treedata_temp(change_node-1,2)=split_var;
  treedata_temp(change_node-1,3)=split_point;
    
//#find left daughter and update tree_matrix and mean st dev if it is a terminal node
  left_daughter=treedata_temp(change_node-1,0);
  right_daughter=treedata_temp(change_node-1,1);
  
  
  IntegerVector nodes_to_update;
  nodes_to_update.push_back(change_node);
 while(ind_nodes_to_update!=0){ 
  for(int i=0;i<ind_nodes_to_update;i++){
    
    NumericVector obs_to_update=find_term_obs(tree_matrix_temp,nodes_to_update[i]);
    
    NumericVector term_cols=unique(find_term_cols(tree_matrix_temp,nodes_to_update[i]));
         
      
    for(int j =0;j<obs_to_update.size();j++){
      tree_matrix_temp(obs_to_update[j],term_cols[0]+1)=999;
    }
      
    
    left_daughter=treedata_temp(nodes_to_update[i]-1,0);    
    right_daughter=left_daughter+1;
    split_var=treedata_temp(nodes_to_update[i]-1,2);
    split_point=treedata_temp(nodes_to_update[i]-1,3);
    ld_obs=get_daughter_obs(xmat,obs_to_update,split_var,split_point)[0];
    rd_obs=get_daughter_obs(xmat,obs_to_update,split_var,split_point)[1];
     
    /*now have the new subsets which go to the left and right daughters of the change node
      need to update these in the tree_matrix
      we can only change internal nodes which means there will always be a child so will
      there will always be at least one column extra in the tree matrix for the children nodes*/
   
    tree_matrix_temp=set_daughter(left_daughter,right_daughter,ld_obs,rd_obs,tree_matrix_temp,term_cols[0]);  
        //now see if ld/rd have children nodes
     if(i==0){
      temp_ind_nodes_to_update=ind_nodes_to_update;
      //*************/
      //need to see what is inseide temp_ind_nodes here!!!
 
      IntegerVector temp_temp;
      temp_nodes_to_update=temp_temp;
      temp_nodes_to_update=nodes_to_update;
 
    }    
    if(treedata_temp(left_daughter-1,4)==-1){
 
    /*when at the terminal node of the branch, change the treedata to update new node positions
    for the swapped tree.*/
      if(ld_obs.size()!=0){  
        
        if(term_cols[0]+2 < tree_matrix_temp.ncol()){
          for(int j=(term_cols[0]+2);j<tree_matrix_temp.ncol();j++){ 
            for(int k=0;k<ld_obs.size();k++){
              tree_matrix_temp(ld_obs[k],j)=0;
            }
          }
        }
        
      }
    }else{
          temp_ind_nodes_to_update=temp_ind_nodes_to_update+1;
          temp_nodes_to_update.push_back(left_daughter);
          //left_daughter2<-treedata_temp[left_daughter,1]
          //print("LD not terminal go again!")
       
    }
    if(treedata_temp(right_daughter-1,4)==-1){
     
          /*when at the terminal node of the branch, change the treedata to update the mean and 
          st dev of tree terminal nodes.*/
          if(rd_obs.size()!=0){
            /*treedata[right_daughter,6]<-mean(y[as.numeric(rownames(rd_obs))])
            treedata[right_daughter,7]<-sd(y[as.numeric(rownames(rd_obs))])
            
            if obs were reassigned and some which were previously non-terminal are now terminal
            then values for these obs in future columns have to be set to zero!*/
            
            if(term_cols[0]+2<tree_matrix_temp.ncol()){
              for(int j=(term_cols[0]+2);j<tree_matrix_temp.ncol();j++){ 
                for(int k=0;k<rd_obs.size();k++){
                  tree_matrix_temp(rd_obs[k],j)=0;
                }
              }
              
            }
          }
        }else{
          temp_ind_nodes_to_update=temp_ind_nodes_to_update+1;
         
          temp_nodes_to_update.push_back(right_daughter); 
          
        }
  }  
  // #need to check this when there is more than one child that needs to be updated
 
     for(int k =0;k<ind_nodes_to_update;k++){
     
        temp_nodes_to_update.erase(k); 
        
      }
 
      IntegerVector temp;
      
      nodes_to_update=temp;
      nodes_to_update=temp_nodes_to_update;
     // temp_ind_nodes_to_update=0;
      temp_nodes_to_update=temp;
      
    
      ind_nodes_to_update=temp_ind_nodes_to_update-ind_nodes_to_update;
  
}    
  ret[0]=treedata_temp;
  ret[1]=tree_matrix_temp;
  return(ret);
}
//########################################################################################//
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

List change_tree(NumericMatrix xmat,IntegerMatrix tree_matrix,NumericMatrix treedata){
  NumericMatrix treedata_temp=clone(treedata);
  IntegerMatrix tree_matrix_temp=clone(tree_matrix);
  List ret(2);
  // List test(2);
  IntegerVector internal_nodes;
  IntegerVector change_node1;
  int change_node;
  IntegerVector split_var1;
  int split_var;
  
  if(treedata_temp.nrow()==1){
   
   ret[0]=treedata_temp;
   ret[1]=tree_matrix_temp;
  
  }else{
    internal_nodes=find_internal_nodes(treedata_temp);
    
    if(internal_nodes.size()>1){
      change_node1=csample_num(internal_nodes,1,false);
      change_node=change_node1[0];
     
    }else{
      change_node=internal_nodes[0];
    }
    //#choose a new split variable and split point at random
    
    IntegerVector xx=seq_len(xmat.ncol());
    
    split_var1=csample_num(xx,1,false);
    split_var=split_var1[0];
    
    //#choose a random split point in the range of the split variable for the observations in 
    //#the chosen change node
    double min_point;
    double max_point;
    min_point=min(xmat(_,split_var-1));
    max_point=max(xmat(_,split_var-1));
   
    NumericVector split_point1=runif(1,min_point,max_point);
    double split_point=split_point1[0];
   
    ret=update_changetree_matrix(xmat,treedata_temp,tree_matrix_temp, split_var,split_point, change_node);
    
  }

  
  //ret[0]=treedata_temp[0];
//  ret[1]=tree_matrix_temp[1];
  return(ret);
}
//############################################################################################//
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector calc_rowsums(NumericMatrix predictions,NumericVector response){
  NumericVector row_sums=response.size();
  NumericVector resids=response.size();
  int ncol=predictions.ncol();
  int nrow=predictions.nrow();
  for(int j=0;j<nrow;j++){
    for(int i=0;i<ncol;i++){  
      row_sums[j]+=predictions(j,i);
    }  
  }
  return(row_sums);
}
//############################################################################################//
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector calculate_resids(NumericMatrix predictions,NumericVector response){
  NumericVector resids=response.size();
  NumericVector row_sums=calc_rowsums(predictions,response);
  resids=response - row_sums;
  return(resids);
}
//#############################################################################################//
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
bool Find_element(int search_node,IntegerVector internal_nodes){
 // if(treedata[internal_nodes[i],1] %in% internal_nodes)
  bool element=std::find(internal_nodes.begin(), internal_nodes.end(),search_node)!=internal_nodes.end();
  return(element);
}
//############################################################################################//

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List append_parent_child(NumericMatrix treedata_temp,IntegerVector internal_nodes){
  List parent_child(2);
    
  //go through each internal node in order find their left and right daughters, if children
  //are also in "internal_nodes" add these to the list of nodes to be swapped
    
  IntegerVector swap_parent;
  IntegerVector swap_child;
 
  for(int i=0; i<internal_nodes.size();i++){
    bool elem_pres=Find_element(treedata_temp(internal_nodes[i]-1,0),internal_nodes);
    
    if(elem_pres){
      swap_parent.push_back(internal_nodes[i]);
      swap_child.push_back(treedata_temp(internal_nodes[i]-1,0));  
    }
     elem_pres=Find_element(treedata_temp(internal_nodes[i]-1,1),internal_nodes);
     if(elem_pres){
        swap_parent.push_back(internal_nodes[i]);
        swap_child.push_back(treedata_temp(internal_nodes[i]-1,1));  
      }
  }
  
  parent_child[0]=swap_parent;
  parent_child[1]=swap_child;
  return(parent_child);
}
//############################################################################################//
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix set_swap(NumericMatrix treedata_temp,int swap_parent,int swap_child){
  double temp_sv;
  double temp_sp;
  temp_sv=treedata_temp(swap_parent-1,2);
  temp_sp=treedata_temp(swap_parent-1,3);
  
  treedata_temp(swap_parent-1,2)=treedata_temp(swap_child-1,2);
  treedata_temp(swap_parent-1,3)=treedata_temp(swap_child-1,3);
  
  treedata_temp(swap_child-1,2)=temp_sv;
  treedata_temp(swap_child-1,3)=temp_sp;
  return(treedata_temp);
}
//############################################################################################//
//***************************************************************************************/

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

IntegerMatrix update_tree_matrix(NumericMatrix xmat,int swap_parent,NumericMatrix treedata_temp,IntegerMatrix tree_matrix_temp){
  int left_daughter;
  int right_daughter;
  left_daughter=treedata_temp(swap_parent-1,0);
  right_daughter=treedata_temp(swap_parent-1,1);
  int ind_nodes_to_update=1;
  IntegerVector nodes_to_update;
  nodes_to_update.push_back(swap_parent);  
  int split_var;
  double split_point;
  IntegerVector ld_obs;
  IntegerVector rd_obs;
  int temp_ind_nodes_to_update=ind_nodes_to_update;
  IntegerVector temp_nodes_to_update;
  
  while(ind_nodes_to_update!=0){
  //go in here when there is still a node that was split so observations in children 
  //nodes need to be reassigned
 
      //go in here when there is still a node that was split so observations in children 
      //nodes need to be reassigned
    for(int i=0;i<ind_nodes_to_update;i++){
 
      NumericVector obs_to_update=find_term_obs(tree_matrix_temp,nodes_to_update[i]);
    
      NumericVector term_cols=unique(find_term_cols(tree_matrix_temp,nodes_to_update[i]));
       
      for(int j =0;j<obs_to_update.size();j++){
        tree_matrix_temp(obs_to_update[j],term_cols[0]+1)=999;
      }
      
      left_daughter=treedata_temp(nodes_to_update[i]-1,0);    
      right_daughter=left_daughter+1;
      split_var=treedata_temp(nodes_to_update[i]-1,2);
      split_point=treedata_temp(nodes_to_update[i]-1,3);
      
      //this means that we haven't got to the end of the branch yet
        //follow it down and keep refining the subset until we get to a terminal node
      
      
      ld_obs=get_daughter_obs(xmat,obs_to_update,split_var,split_point)[0];
      rd_obs=get_daughter_obs(xmat,obs_to_update,split_var,split_point)[1];
     
     /*now have the new subsets which go to the left and right daughters of the change node
        need to update these in the tree_matrix
        we can only change internal nodes which means there will always be a child so will
        there will always be at least one column extra in the tree matrix for the children nodes*/
     
      tree_matrix_temp=set_daughter(left_daughter,right_daughter,ld_obs,rd_obs,tree_matrix_temp,term_cols[0]);  
      
      //if i>1 this is overwriting the nodes to update from previous internal nodes
        //at the same depth, try this
     // if(i==0){
      //  int temp_ind_nodes_to_update=ind_nodes_to_update;
      //  IntegerVector temp_nodes_to_update=nodes_to_update;
    //  }
     //now see if ld/rd have children nodes
        
      //if i>1 this is overwriting the nodes to update from previous internal nodes
        //at the same depth, try this
    if(i==0){
      temp_ind_nodes_to_update=ind_nodes_to_update;
      //*************/
      //need to see what is inseide temp_ind_nodes here!!!
  
      IntegerVector temp_temp;
      temp_nodes_to_update=temp_temp;
      temp_nodes_to_update=nodes_to_update;
     
    }

    //now see if ld/rd have children nodes
    if(treedata_temp(left_daughter-1,4)==-1){

    /*when at the terminal node of the branch, change the treedata to update new node positions
    for the swapped tree.*/
      if(ld_obs.size()!=0){  
        
        if(term_cols[0]+2 < tree_matrix_temp.ncol()){
          for(int j=(term_cols[0]+2);j<tree_matrix_temp.ncol();j++){ 
            for(int k=0;k<ld_obs.size();k++){
              tree_matrix_temp(ld_obs[k],j)=0;
            }
          }
        }
        
      }
    }else{
          temp_ind_nodes_to_update=temp_ind_nodes_to_update+1;
          temp_nodes_to_update.push_back(left_daughter);
       
    } 
    if(treedata_temp(right_daughter-1,4)==-1){
     
          /*when at the terminal node of the branch, change the treedata to update the mean and 
          st dev of tree terminal nodes.*/
          if(rd_obs.size()!=0){
            /*treedata[right_daughter,6]<-mean(y[as.numeric(rownames(rd_obs))])
            treedata[right_daughter,7]<-sd(y[as.numeric(rownames(rd_obs))])
            
            if obs were reassigned and some which were previously non-terminal are now terminal
            then values for these obs in future columns have to be set to zero!*/
            
            if(term_cols[0]+2<tree_matrix_temp.ncol()){
              for(int j=(term_cols[0]+2);j<tree_matrix_temp.ncol();j++){ 
                for(int k=0;k<rd_obs.size();k++){
                  tree_matrix_temp(rd_obs[k],j)=0;
                }
              }
              
            }
          }
        }else{
          temp_ind_nodes_to_update=temp_ind_nodes_to_update+1;
        
          temp_nodes_to_update.push_back(right_daughter); 
        
          
        }
  }
  // #need to check this when there is more than one child that needs to be updated

     for(int k =0;k<ind_nodes_to_update;k++){

        temp_nodes_to_update.erase(k); 
        
      }

      IntegerVector temp;
      
      nodes_to_update=temp;
      nodes_to_update=temp_nodes_to_update;
     // temp_ind_nodes_to_update=0;
      temp_nodes_to_update=temp;
      

      ind_nodes_to_update=temp_ind_nodes_to_update-ind_nodes_to_update;
      
}    
             
  return(tree_matrix_temp);
}
//***************************************************************************************//

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

List swap_nodes(NumericMatrix xmat,NumericMatrix treedata,IntegerMatrix tree_matrix){
  IntegerVector internal_nodes;
  NumericMatrix treedata_temp=clone(treedata);
  IntegerMatrix tree_matrix_temp=clone(tree_matrix);
  internal_nodes=find_internal_nodes(treedata_temp);
  
  List parent_child=append_parent_child(treedata_temp,internal_nodes);
  
  IntegerVector parent=parent_child[0];
  IntegerVector child=parent_child[1];
  double yy=parent.size();
  IntegerVector xx=seq_len(yy);
  int swap_parent;
  int swap_child;
  List swap_output(2);
  IntegerMatrix tree_matrix2;
  NumericMatrix treedata2;
  
 
  if((treedata_temp.nrow()==1) | (tree_matrix_temp.ncol()==2)){
    swap_output[0]=treedata_temp;
    swap_output[1]=tree_matrix_temp;
  }else{
    IntegerVector swap_node=csample_num(xx,1,false);
    int swap_node1=swap_node[0];
    swap_parent=parent[swap_node1-1];
    swap_child=child[swap_node1-1];
    
   //create set_swap() function
    treedata2=set_swap(treedata_temp,swap_parent,swap_child); 
   
    //find left daughter and update tree_matrix and mean st dev if it is a terminal node
    tree_matrix2=update_tree_matrix(xmat,swap_parent,treedata_temp,tree_matrix_temp);
    
    swap_output[0]=treedata_temp;
    swap_output[1]=tree_matrix_temp;
  }
  
 return(swap_output);
}

//############################################################################################//
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
//###########################################################################################//
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
//###########################################################################################//
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
  IntegerVector find_term_child(NumericMatrix tree_table){
  NumericVector node_status=tree_table.column(4);
  IntegerVector internal_nodes_prop;
  
  for(int i=0;i<node_status.size();i++){
  
    if(node_status[i]==1){
  
      internal_nodes_prop.push_back(i+1);
      
    }
  }
  
  IntegerVector both_child_term(internal_nodes_prop.size());
 
  for(int k=0;k<internal_nodes_prop.size();k++){
 
    if((tree_table(tree_table(internal_nodes_prop[k]-1,0)-1,4)==
    tree_table(tree_table(internal_nodes_prop[k]-1,1)-1,4)) && (tree_table(tree_table(internal_nodes_prop[k]-1,0),4)==-1)){
 
      both_child_term[k]=1;
    
    }else{
      both_child_term[k]=0;
    }
  }

  IntegerVector nodes_to_prune;
  
  for(int k=0;k<both_child_term.size();k++){
    if(both_child_term[k]==1){
      nodes_to_prune.push_back(internal_nodes_prop[k]);     
    }
  }  

  return(nodes_to_prune);
}
    
//#######################################################################################  
#include <Rcpp.h>
//#include "find_term_child.h"
using namespace Rcpp;
// [[Rcpp::export]]
double choose_prune(double ratio,double tree_ratio,NumericMatrix prior_tree,NumericMatrix proposal_tree){
  double ratio_out;
  double var_density=1;
  double q_t_tstar;
  double q_tstar_t;
  double qratio;
  
  NumericVector terminal_nodes_pro;
  IntegerVector nodes_to_prune=find_term_child(prior_tree);

  if(prior_tree.nrow()==1){
    q_t_tstar=1;

  }else{
    double node_size=nodes_to_prune.size();
    double num_nodes=1/node_size;
    q_t_tstar=0.25*num_nodes;

  }
  
  for(int l=0;l<proposal_tree.nrow();l++){
    
    if(proposal_tree(l,4)==-1){
      
      terminal_nodes_pro.push_back(l+1);
    }
  }
  //for(int i=0;i<terminal_nodes_pro.size();i++)
  q_tstar_t=0.25*(1.0/terminal_nodes_pro.size())*(var_density);
  qratio=q_tstar_t/q_t_tstar;
  ratio_out=ratio*tree_ratio*qratio;

  return(ratio_out);
}
//##########################################################################################//
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
double choose_grow(double ratio,double tree_ratio,NumericMatrix prior_tree,NumericMatrix proposal_tree){
  double var_density=1;
  double ratio_out;
  double q_t_tstar;
  double q_tstar_t;
  double qratio;
  NumericVector terminal_nodes_pri;
  IntegerVector nodes_to_prune;
  
  for(int l=0;l<prior_tree.nrow();l++){
    
    if(prior_tree(l,4)==-1){
      
      terminal_nodes_pri.push_back(l+1);
      
    }
  }
  double term_node_pri=terminal_nodes_pri.size();
 
  q_t_tstar=0.25*(1/term_node_pri)*var_density;
  //get probability of going from proposal tree to prior tree which is
  //p(prune)*p(prune node chosen)
  
  nodes_to_prune=find_term_child(proposal_tree);
  
  //p(prune)*p(prune node chosen)
  double prune_nodes=nodes_to_prune.size();
  q_tstar_t=0.25*(1/prune_nodes);
  
  qratio=q_tstar_t/q_t_tstar; 

  ratio_out=ratio*tree_ratio*qratio;

  return(ratio_out);
}
//#########################################################################################//
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
//#########################################################################################//
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector get_new_mean(NumericVector terminal_nodes,List new_mean_var){
  NumericVector node_means;
  for(int k=0;k<terminal_nodes.size();k++){
    NumericVector sd=new_mean_var[1];
    NumericVector temp_mean=new_mean_var[0];
    double new_mean= R::rnorm(temp_mean[k],sqrt(sd[k]));
    node_means.push_back(new_mean);
   
  }

  return(node_means);
}
//##########################################################################################//      
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
double update_sigma(double a1,double b,NumericVector resids,int n){
  NumericVector sq_resid=resids*resids;
  double ssr=0;
  for(int i=0;i<resids.size();i++){
    ssr+=sq_resid[i];
  }

  double shape=(a1+n/2);
  double rate =((ssr/2)+(1/b));

  RNGScope scope;
  double tau =R::rgamma(shape,1/rate);
  double sigma = sqrt((1/tau));
  return sigma;

}
//######################################################################################//
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector get_original(double low,double high,double sp_low,double sp_high,NumericVector sum_preds){
  NumericVector original_y=(sum_preds*(-low+high))/(-sp_low+sp_high) + (-high*sp_low+low*sp_high)/(-sp_low+sp_high);
  return(original_y);
}
//######################################################################################//
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector get_imp_vars(NumericVector tree_table,int num_col,NumericVector current_vars){
  
  NumericVector vars_chosen=sort_unique(tree_table);
  
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
//#########################################################################################//
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

List BART_Main(NumericMatrix data,NumericVector y,int niter,int num_trees,double sigma,double nu,double lambda,double a,double alpha,double beta){
  int n=data.nrow();
 // int vardensity=1;
  bool node_grown;
  bool node_pruned;
  NumericVector y_scaled=scale_response(min(y),max(y),-0.5,0.5,y=y);

  double mu_mu=mean(y_scaled)/num_trees;
  double sigma_mu=sd(y_scaled)/num_trees;
  
   //initialize trees
  
  List prior_tree=start_tree(num_trees,mu_mu,sigma_mu);
  List prior_matrix=start_matrix(num_trees,n);

  List post_tree=prior_tree;
  List post_matrix=prior_matrix;
  
  //  #make start predictions
  NumericMatrix predy(n,num_trees);
  predy=start_predy(num_trees,n=n,mu_mu,sigma_mu);
  
  NumericVector sigma_its(niter);
  
  double sigma_mu_sq=pow(sigma_mu,2);
  
  // #initialize variables
  NumericMatrix var_imp(niter,data.ncol());
  NumericMatrix post_predictions(niter,y_scaled.size());
  NumericMatrix post_predictions_orig=post_predictions;
  NumericVector sum_preds;
  NumericVector original_y;
  NumericVector pri(num_trees);
  NumericMatrix proposal_table;
  IntegerMatrix proposal_matrix;
  List mu_term_nodes(num_trees);
  double prop_lik;
  double ratio;
  double tree_ratio=0;
  double a1;
  double b;
  NumericVector sd;
  List proposal_tree;
  NumericVector terminal_nodes;
  List updated_preds;
 
  NumericVector temp_mean;
  List new_mean_var;
  
  for(int j=0;j<niter;j++){

    for(int i=0;i<num_trees;i++){
       node_pruned=0;
       node_grown=0;
     //get predicted y for all except tree i
        NumericMatrix test_predy=clone(predy);
        NumericMatrix temp_predy=remove_curr_col(test_predy,i);
       
        //create a temperary matrix and delete current col (change to armadillo mat)
        NumericVector resids=calculate_resids(temp_predy,y_scaled);
      
        //get likelihood of the current tree
        NumericMatrix temp_tree=post_tree[i];
        IntegerMatrix temp_mat=post_matrix[i];
      
        pri[i]=likelihood_function(resids,temp_tree,temp_mat,a,mu_mu,nu,lambda);
        
        //now get proposal tree  
        //generate random number and choose between change and swap steps!
        NumericVector rand=runif(1);
      
        if(rand[0]<=0.25){
          //choose change step
      
          NumericMatrix test_tree=prior_tree[i];
          IntegerMatrix test_mat=prior_matrix[i];
          proposal_tree=change_tree(data,prior_matrix[i],prior_tree[i]);

        }else if((rand[0]>0.25) & (rand[0]<=0.5)){
          //choose swap step
              NumericMatrix test_tree=prior_tree[i];
          IntegerMatrix test_mat=prior_matrix[i];
          proposal_tree=swap_nodes(data,prior_tree[i],prior_matrix[i]);
        
        }else if((rand[0]>0.5) & (rand[0]<=0.75)){
          //choose grow
          IntegerVector terminal_node;
          NumericMatrix tree_table=prior_tree[i];
           for(int k=0;k<tree_table.nrow();k++){
              if(tree_table(k,4)==-1){
                terminal_node.push_back(k+1);
              }
            }
            
          NumericMatrix test_tree=prior_tree[i];
          IntegerMatrix test_mat=prior_matrix[i];
          IntegerVector grow_node=csample_num(terminal_node,1,false);
          proposal_tree=grow_tree(data,resids,prior_matrix[i],grow_node[0],prior_tree[i]);
         
          node_grown=1;
        }else{
          //choose prune
         
          node_pruned=1;
          proposal_tree=prune_tree(data,prior_matrix[i],prior_tree[i]);
        }
  
      prop_lik=likelihood_function(resids,proposal_tree[0],proposal_tree[1],a,mu_mu,nu,lambda);
      ratio=exp(prop_lik-pri[i]);

      if((node_grown==1) | (node_pruned==1)){
       
        tree_ratio=choose_grow_prune(proposal_tree[0],proposal_tree[1],prior_tree[i],prior_matrix[i],alpha,beta);

      }
    
    if(node_grown==1){
      //get q(T'-> T)
      //get prob of going from prior tree to proposal tree which is
      //p(grow)*p(grow node chosen)*p(split var chosen)*p(split point chosen)
      
      ratio=choose_grow(ratio,tree_ratio,prior_tree[i],proposal_tree[0]);
      
    }
    if(node_pruned==1){
      //get probability of going from proposal tree to prior tree which is
      //p(prune)*p(prune node chosen)    
      ratio=choose_prune(ratio,tree_ratio,prior_tree[i],proposal_tree[0]);
    }
    //now choose which tree to update
    //now update step with the most likely tree
    
    if(ratio<1){
//      #if ratio<1 then set tree as proposal tree with probability ratio else keep the prior tree
  NumericVector    prob=runif(1);
//      #probs[i]<-prob
      if(prob[0]<ratio){
      // #if random number is less than ratio use proposal tree
        post_tree[i]=proposal_tree[0];
        prior_tree[i]=post_tree[i];
        post_matrix[i]=proposal_tree[1];
        prior_matrix[i]=post_matrix[i];
        //#here I need to save the yhat and importance variables for each tree
        //# terminal_nodes<-which(list_of_tree_tables[[i]][[j]][,5]==-1)
        terminal_nodes=find_term_nodes(post_tree[i]);        
      }else{
        //#else use prior
        
        post_tree[i]= prior_tree[i];
        
        post_matrix[i]=prior_matrix[i];
      }
    }else{
      //#if ratio>1 use proposal
      
      post_tree[i]=proposal_tree[0];
      post_matrix[i]=proposal_tree[1];
      prior_tree[i]=post_tree[i];
      terminal_nodes=find_term_nodes(post_tree[i]);

      prior_matrix[i]=post_matrix[i];
      
    }
      //terminal_nodes<-which(list_of_tree_tables[[i]][[j]][,5]==-1)
      terminal_nodes=find_term_nodes(post_tree[i]); 
   
    new_mean_var=update_mean_var(post_tree[i],post_matrix[i],resids,a,sigma,mu_mu);
   NumericVector new_mean=get_new_mean(terminal_nodes,new_mean_var);
     
   updated_preds=update_predictions(post_tree[i],post_matrix[i],
                                      new_mean,new_mean_var[1],n);
    
    NumericVector temp_preds=updated_preds[1];
       

        predy(_,i)=temp_preds;  
    
      post_tree[i]=updated_preds[0];
 
//      #now have chosen tree add the variable counts
//    #vars chosen in posterior tree
//    # var_imp[j,]<-get_imp_vars(tree_table=list_of_tree_tables[[i]][[j]][,3],num_col=ncol(data),current_vars=var_imp[j,])
      NumericMatrix tree=post_tree[i];
      NumericVector impvar=tree(_,2);
      var_imp(j,_)=get_imp_vars(impvar,data.ncol(),var_imp(j,_));
    
    
  //  #update tau for next iteration
  //  #set a and b (changing parameters from (nu.lambda) in inverse gamma to gamma)
    a1=nu/2;
    b=2/(nu*lambda);
    
//    #calculate residuals (response for tree)
//    #sum_preds<-rep(0,n)
    
  //  #sum_preds<-apply(predy,1,sum)
    resids=calculate_resids(predy,y_scaled);
    sum_preds=calc_rowsums(predy,y_scaled);

      sigma=update_sigma(a1,b,resids,n);   

//    #update residuals for next iteration
    sigma_its[j]=sigma; 
       
    }
    post_predictions(j,_)=calc_rowsums(predy,y_scaled);
    sum_preds=calc_rowsums(predy,y_scaled);
    original_y=get_original(min(y),max(y),-0.5,0.5,sum_preds);
    
    post_predictions_orig(j,_)=original_y;
  }
  List ret(5);
  ret[0]=post_tree;
  ret[1]=post_matrix;
  ret[2]=sigma_its;
  ret[3]= var_imp;
  ret[4]=post_predictions_orig;
  return(ret);
}