
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