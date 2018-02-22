library(Rcpp)

cppFunction('
  int add(int x, int y, int z) {
    int sum = x + y + z;
    return sum;
  }'
)
add # like a regular R function, printing displays info about the function
#> function (x, y, z) 
#> .Primitive(".Call")(<pointer: 0x107902ba0>, x, y, z)
#> <environment: 0x1034b2268>
test<-add(1, 2, 30)
test
#> [1] 6

########################################################
cppFunction('
  int one() {
    return 1;
  }
')
one()
########################################################
signR <- function(x) {
  if (x > 0) {
    1
  } else if (x == 0) {
    0
  } else {
    -1
  }
}

cppFunction('
  int signT(int x){
    if(x>0){
      return 1;
    }else if(x==0){
      return 0;
    }else{
      return -1;
    }
  }
')
signT(-4)
#####################################################################################
sumR <- function(x) {
  total <- 0
  for (i in seq_along(x)) {
    total <- total + x[i]
  }
  total
}

cppFunction('
  double sumC(NumericVector x){
    double total =0;
    int n=x.size();
    for(int i=0;i<n;i++){
      total +=x[i];
    }
    return total;
  }
')

library(microbenchmark)
x <- runif(1e3)
microbenchmark(
  sum(x),
  sumR(x),
  sumC(x)
)
sumC(c(1,2,3))
#######################################################################
pdistR <- function(x, ys) {
  sqrt( (x - ys) ^ 2 )
}
pdistR(c(2),c(3,2,1))

cppFunction('
  NumericVector pdistC(double x,NumericVector ys){
  int n = ys.size();
  NumericVector answer(n);
    for(int i=0;i<n;i++){
      answer[i]=sqrt(pow((x-ys[i]),2));
    }
  return answer;
  }
            ')
pdistC(c(3.2),c(3,2,1))
#######################################################################
#get row sums of a matrix and output vector of sum

cppFunction('
  NumericVector rowsumsC(NumericMatrix x){
    int ncolx=x.ncol();
    std::cout<<"ncolx is"<<ncolx;
    int nrowx=x.nrow();
std::cout<<"nrowx is"<<nrowx;
    NumericVector out(nrowx);
std::cout<<"out is"<<out[1];
    NumericVector rowsums(ncolx);
    for(int i=0;i<nrowx;i++){
        double sumrow=0;
        for(int j=0;j<ncolx;j++){
            sumrow+=x(i,j);
        }
      out[i]=sumrow;
    }
    return(out);
  }
            ')

x <- matrix(sample(100), 10)
rowSums(x)
rowsumsC(x)
###################################################################################
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double meanC(NumericVector x) {
  int n = x.size();
  double total = 0;
  
  for(int i = 0; i < n; ++i) {
    total += x[i] / n;
  }
  return total;
}

sourceCpp("C:/Users/Belinda/Desktop/second_year/BART/BARTC/test_sourcecpp.cpp")
x <- runif(1e5)
microbenchmark(
  mean(x),
  meanC(x))

sourceCpp("C:/Users/Belinda/Desktop/second_year/BART/BARTC/test_pdist_sugar.cpp")
pdistC2(2,c(1,2,3))

####################################################################################
test1<-sourceCpp("C:/Users/Belinda/Desktop/second_year/BART/BARTC/missing_calues_c.cpp")
test1[[1]]

####################################################################################
#                     BART Functions in C++      
####################################################################################

sourceCpp("C:/Users/Belinda/Desktop/second_year/BART/BARTC/predy.cpp")
predy2<-start_predy(num_trees=num_trees,n=10,start_mean=0,start_sd=100)

###############################################################################

sourceCpp("C:/Users/Belinda/Desktop/second_year/BART/BARTC/y_scaled.cpp")
y_scaled2<-scale_response(a=min(y),b=max(y),c=-0.5,d=0.5,y=y)

###############################################################################
sourceCpp("C:/Users/Belinda/Desktop/second_year/BART/BARTC/get_original.cpp")
original_y3<-get_original(low=min(y),high=max(y),sp_low=-0.5,sp_high=0.5,sum_preds=sum_preds)
