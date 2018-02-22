library(Rcpp)
treetable=post_tree[[3]]
treemat=post_matrix[[3]]

sourceCpp("C:/Users/Belinda/Desktop/second_year/BART/BARTC/get_best_likelihood.cpp")
start=Sys.time()
test=get_best_split(resids=resids,data=as.matrix(data),treetable=treetable,treemat=treemat,a=a,mu=mu_mu,nu=nu, lambda=lambda,grid_length=50)
end=Sys.time()
end-start
test[1]
test[2]
test[3]


