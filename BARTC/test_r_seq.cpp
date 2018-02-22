#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

NumericVector create_cut_points(NumericVector x,int lengthout) {
   double maxpt=max(x);
   double minpt=min(x);
   std::cout<<maxpt<<"\n";
   double rangex=maxpt-minpt;
   double increment=rangex/(lengthout-1);
   NumericVector ret(lengthout);
   ret[0]=minpt;
   ret[49]=maxpt;
   for(int i=1;i<49;i++){
    ret[i]=ret[i-1]+increment;        
   }
   return ret;
}
