#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
double test_power_fn(double prior_p,double beta){
  double propsplit=pow(prior_p,-beta);
  return(propsplit);
}