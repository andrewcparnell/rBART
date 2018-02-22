#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
double update_sigma(double a1,double b,NumericVector resids,int n){
  NumericVector sq_resid=resids*resids;
  double ssr=0;
  for(int i=0;i<resids.size();i++){
    ssr+=sq_resid[i];
  }
 // std::cout<<ssr;
  double shape=(a1+n/2);
  double rate =((ssr/2)+(1/b));
  //std::cout<<"shape"<<shape<<" "<<"rate"<<rate;
  RNGScope scope;
  double tau =R::rgamma(shape,1/rate);
  double sigma = sqrt((1/tau));
  return sigma;
//  return sigma;
}