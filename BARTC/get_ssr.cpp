#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
double get_ssr(NumericVector resids){
  NumericVector sq_resid=resids*resids;
  double ssr=0;
  for(int i=0;i<resids.size();i++){
    ssr+=sq_resid[i];
  }
  return(ssr);
}
/*
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
double gamma_rand(){
  //double rate=b*3;
  double tau =R::rgamma(5,1/4);
  return(tau);
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericVector find_term_child(NumericMatrix tree_table){
  NumericVector node_status=tree_table.column(4);
  IntegerVector internal_nodes_prop=std::which(node_status==1);
  return(node_status);
}*/