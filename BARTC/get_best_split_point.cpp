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
  std::cout<<"x is"<<x.size()<<" "<<x[0]<<"\n";
  IntegerVector ret = RcppArmadillo::sample(x, size, replace, prob);
  return ret;
}
//######################################################################################//
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
//####################################################################################//

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

NumericVector create_cut_points(NumericVector x,int lengthout) {
   double maxpt=max(x);
   double minpt=min(x);
   double rangex=maxpt-minpt;
   double increment=rangex/(lengthout-1);
   NumericVector ret(lengthout);
   ret[0]=minpt;
   ret[lengthout-1]=maxpt;
   for(int i=1;i<(lengthout-1);i++){
    ret[i]=ret[i-1]+increment;        
   }
   return ret;
}
//####################################################################################//
NumericVector get_best_split(NumericVector resids,NumericMatrix treetable,IntegerMatrix treemat,double a,double mu,double nu,double lambda,int grid_length){
  //this function will calculate the highest likelihood for different values of
  //split points of x
  NumericMatrix treetable_c=clone(treetable);
  IntegerMatrix treemat_c=clone(treemat);
  //find the internal nodes in the tree
  IntegerVector split_vars=find_internal_nodes(treetable_c);
  
  //randomly choose a split variable from list of internal nodes
  IntegerVector sp=csample_num(split_vars,1,false);
  
  //split x into 50 grid points  
  NumericVector grid_points=create_cut_points(sp[0],grid_length);
  //for each grid point have to change the split points of the tree and reassign subjects to terminal nodes
  //in order to calculate the likelihood for the slightly changed tree
  
  
  //lik=likelihood_function(resids,change_treetable,change_treemat,a,mu_mu,nu,lambda);
  
  return grid_points;
}