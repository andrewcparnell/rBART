#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp ;

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
    std::cout<<sv_col[obs_to_update[i]-1]<<"\n";
    
    if(sv_col[obs_to_update[i]-1]<=split_point){
      std::cout<<"xmat col2"<<sv_col[obs_to_update[i]-1]<<"\n";
     
      ld_obs.push_back(obs_to_update[i]-1);
      //std::cout<<"ld_obs"<<ld_obs[0]<<"\n";
    }
    if(sv_col[obs_to_update[i]-1]>=split_point){
      std::cout<<"xmat col2"<<sv_col[obs_to_update[i]-1]<<"\n";
     
      rd_obs.push_back(obs_to_update[i]-1);
      //std::cout<<"ld_obs"<<ld_obs[0]<<"\n";
    }
  }
  //ld_obs=xmat(obs_to_update[,1],)
  //ld_obs<-ld_obs[which(ld_obs[,split_var]<=split_point),]
  daughter_obs[0]=ld_obs;
  daughter_obs[1]=rd_obs;
  return(daughter_obs);
}