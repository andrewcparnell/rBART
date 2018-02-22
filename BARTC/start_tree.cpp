#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List start_tree(int num_trees,double start_mean,double start_sd){
  List pri_tree(num_trees);
  NumericMatrix treemat(1,7);
  
  
  for(int i =0;i<num_trees;i++){
    double rand=R::rnorm(start_mean,start_sd);
    //std::cout<<rand;
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