
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
//******************************************************************************************//
//get the log likelihood
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

IntegerVector likelihood_fun(NumericVector y,NumericMatrix treetable,IntegerMatrix obs_to_nodes,double a,double mu,double nu,double lambda){
  
  IntegerVector terminal_nodes;

  for(int l=0;l<treetable.nrow();l++){
    
    if(treetable(l,4)==-1){
      terminal_nodes.push_back(l+1);
   //   std::cout<<"i is"<<l<<"\n";      
    }
  }
  int b=terminal_nodes.size();
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
  
  for(int i=0;i< b;i++){
    //number of observations in terminal nodes
    term_obs=find_term_obs(obs_to_nodes,terminal_nodes[i]);
    ni=term_obs.size();
    n[i]=ni;
    
    std::cout<<"ni is"<<n[i]<<"\n";
    double temp1=0;
    double temp2=0;
   
    
    for(int j=0;j<ni;j++){
      temp1+=y[term_obs[j]];    
      temp2+=pow(y[term_obs[j]],2);    
    }
     std::cout<<"temp1"<<temp1<<"\n";
  std::cout<<"temp2"<<temp2<<"\n";
   sumr[i]=temp1;
    sumsqr[i]=temp2;
    s[i] = (sumsqr[i]-pow(sumr[i],2)/n[i]);
     //s=s[!is.na(s)]
    t[i] = (n[i]*a/(n[i]+a))*pow((sumr[i]/n[i]-mu),2);
     //t=t[!is.na(t)]
     std::cout<<"s"<<t[i]<<"\n";
  std::cout<<"t"<<s[i]<<"\n";
  } 
  double p1=(b/2)*log(a);
  double p2=sum(0.5*log(n+a));
  double p3=(sum(n)+nu)/2*log(sum((s+t)+nu*lambda));
  double tree_log_lik=p1 - p2 -p3;
  //double tree_log_lik=(b/2)*log(a) - sum(0.5*log(n+a)) -(sum(n)+nu)/2*log(sum((s+t)+nu*lambda));
 std::cout<<"p1"<<p1<<"\n"; 
 std::cout<<"p2"<<p2<<"\n"; 
 std::cout<<"p3"<<p3<<"\n"; 
  std::cout<<"log(5)"<<log(5)<<"\n"; 
 std::cout<<"tree log lik"<<tree_log_lik<<"\n"; 
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
return(terminal_nodes);
}