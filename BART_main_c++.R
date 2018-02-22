####################################################################################
####################################################################################
library(Rcpp)

sourceCpp("BARTC/predy.cpp")
sourceCpp("BARTC/y_scaled.cpp")
sourceCpp("BARTC/get_original.cpp")
sourceCpp("BARTC/calculate_resids.cpp")
sourceCpp("BARTC/start_matrix.cpp")
sourceCpp("BARTC/start_tree.cpp")
sourceCpp("BARTC/update_sigma.cpp")
##sourceCpp("BARTC/find_term_child.cpp")
sourceCpp("BARTC/choose_prune.cpp")
sourceCpp("BARTC/choose_grow.cpp")
##sourceCpp("BARTC/get_tree_prior.cpp")
sourceCpp("BARTC/choose_grow_prune.cpp")
sourceCpp("BARTC/update_mean_var.cpp")
sourceCpp("BARTC/update_predictions.cpp")
sourceCpp("BARTC/get_imp_vars.cpp")
sourceCpp("BARTC/likelihood_function.cpp")
sourceCpp("BARTC/swap_nodes.cpp")
sourceCpp("BARTC/grow_tree.cpp")
sourceCpp("BARTC/prune_tree.cpp")
sourceCpp("BARTC/change_nodes.cpp")

set.seed(11)
num_trees<-3
niter<-5000
#set priors muij and sigmamu such that p(muij|Tj)~N(num_trees*mu,num_trees*sd(y))
#as per BART paper here num_trees=2
n=nrow(data)
#first scale y 
nu=3
lambda=0.1
y_scaled<-scale_response(a=min(y),b=max(y),c=-0.5,d=0.5,y=y)

mu_mu<-mean(y_scaled)/num_trees
sigma_mu<-sd(y_scaled)/num_trees

#write function to initialize trees
#list_of_tree_tables<-start_tree1(y=y_scaled,num_trees=num_trees,start_mean=mu_mu,start_sd=sigma_mu)
#list_of_tree_matrix<-start_matrix1(y=y_scaled,num_trees=num_trees,n=n)
prior_tree<-start_tree(num_trees=num_trees,start_mean=mu_mu,start_sd=sigma_mu)
prior_matrix<-start_matrix(num_trees=num_trees,n=n)
#for(i in 1:num_trees){
#  prior_tree[[i]]<-list_of_tree_tables[[i]][[1]]
#  prior_matrix[[i]]<-list_of_tree_matrix[[i]][[1]]
post_tree<-prior_tree
post_matrix<-prior_matrix

#############################################################################################
start=Sys.time() 
#set.seed(101)
#make start predictions
predy<-start_predy(num_trees=num_trees,n=n,start_mean=mu_mu,start_sd=sigma_mu)

sigma_its=rep(0,niter)

sigma_mu_sq<-sigma_mu^2

#parameters to choose
#set global parameters

a=1/3

#alpha and beta are shape and depth parameters for tree prior
alpha=0.95
beta=2
var_density<-1

#initialize variables
var_imp<-matrix(0,nrow=niter,ncol=ncol(data))
post_predictions<-matrix(nrow=niter,ncol=length(y_scaled))
post_predictions_orig<-post_predictions
post_predictions2<-matrix(nrow=niter,ncol=length(y_scaled))
post_predictions_orig2<-post_predictions

pri<-rep(0,num_trees)

pri2<-rep(0,num_trees)
mu_term_nodes<-list()
length(mu_term_nodes)<-num_trees
#sigma<-sigma_true
#sigma<-rgamma(1,shape=1,rate=1)
templm = lm(y_scaled~.,data=as.data.frame(data)) 
sigma = summary(templm)$sigma
#sigma=sigma_true
for(j in 2:niter){
  if(j%%1000==0) print(c("Number of iterations ",j))
  for(i in 1:num_trees){
    node_pruned=node_grown=0
    #get predicted y for all except tree i
    #resids_s<-calculate_resids(predictions=predy[,-i],response=y_scaled) 
    resids<-calculate_resids(predictions=predy[,-i],response=y_scaled) 
    #resids<-resids_s[[1]]
    #list_of_tree_tables[[i]][[j]]<-list_of_tree_tables[[i]][[j-1]]
    #post_tree[[i]]<-prior_tree[[i]]
    #list_of_tree_matrix[[i]][[j]]<-list_of_tree_matrix[[i]][[j-1]]
    #post_matrix[[i]]<-prior_matrix[[i]]
    #get likelihood of the current tree
    #pri[i]<-likelihood_function(xmat=data,y=resids,treetable=list_of_tree_tables[[i]][[j-1]],obs_to_nodes=list_of_tree_matrix[[i]][[j-1]],nu,lambda,mu=mu_mu,a=1/3)
    pri[i]<-likelihood_function(y_temp=resids,treetable_temp=post_tree[[i]],obs_to_nodes_temp=post_matrix[[i]],nu,lambda,mu=mu_mu,a=1/3)
    #now get proposal tree  
    #generate random number and choose between change and swap steps!
    rand<-runif(1)
    rand
    if(rand<=0.25){
      empty_nodes=0
      
      #proposal_tree<-change_tree(xmat=data,y=resids,empty_nodes=empty_nodes,tree_matrix=list_of_tree_matrix[[i]][[j-1]],treedata=list_of_tree_tables[[i]][[j-1]])
      proposal_tree<-change_tree(xmat=as.matrix(data),tree_matrix=prior_matrix[[i]],treedata=prior_tree[[i]])
      proposal_table=proposal_tree[[1]]
      proposal_matrix=proposal_tree[[2]]
    }else if(rand>0.25 & rand<=0.5){
      
      #proposal_tree<-swap_nodes(xmat=data,y=resids,tree_matrix=list_of_tree_matrix[[i]][[j-1]],treedata=list_of_tree_tables[[i]][[j-1]])
      proposal_tree<-swap_nodes(xmat=as.matrix(data),tree_matrix=prior_matrix[[i]],treedata=prior_tree[[i]])
      proposal_table=proposal_tree[[1]]
      proposal_matrix=proposal_tree[[2]]
    }else if(rand>0.5 & rand<=0.75){
      
      #      terminal_node<-which(list_of_tree_tables[[i]][[j-1]][,5]==-1)
      terminal_node<-which(post_tree[[i]][,5]==-1)
      grow_node<-terminal_node[sample(1:length(terminal_node),1)]
      #  set.seed(100)
      # proposal_tree<-grow_tree(xmat=data,y=resids,prior_tree_matrix=list_of_tree_matrix[[i]][[j-1]],grow_node=grow_node,prior_tree_table=list_of_tree_tables[[i]][[j-1]])
      proposal_tree<-grow_tree(xmat=as.matrix(data),y=resids,prior_tree_matrix=prior_matrix[[i]],grow_node=grow_node,prior_tree_table=prior_tree[[i]])
      proposal_table=proposal_tree[[1]]
      proposal_matrix=proposal_tree[[2]]
      node_grown<-1
    }else{
      
      #proposal_tree<-prune_tree(xmat=data,y=resids,prior_tree_matrix=list_of_tree_matrix[[i]][[j-1]],prior_tree_table=list_of_tree_tables[[i]][[j-1]])      
    #  proposal_tree<-prune_tree(xmat=data,y=resids,prior_tree_matrix=prior_matrix[[i]],prior_tree_table=prior_tree[[i]])      
      node_pruned<-1
      proposal_tree<-prune_tree(xmat=as.matrix(data),prior_tree_matrix=prior_matrix[[i]],prior_tree_table=prior_tree[[i]])
      proposal_table=proposal_tree[[1]]
      proposal_matrix=proposal_tree[[2]]
    }
    
    prop_lik<-likelihood_function(y_temp=resids,treetable_temp=proposal_table,obs_to_nodes_temp=proposal_matrix,nu=nu,lambda=lambda,mu=mu_mu,a=1/3)
    
    # pro[[i]]<-prop_lik
    #ratio<-exp(prop_lik-pri[[i]])
    ratio<-exp(prop_lik-pri[[i]])
    
    if(node_grown==1 | node_pruned==1){
      #tree_ratio<-choose_grow_prune(proposal_table=proposal_table,proposal_matrix=proposal_matrix,prior_table=list_of_tree_tables[[i]][[j-1]],prior_matrix=list_of_tree_matrix[[i]][[j-1]],alpha=alpha,beta=beta)
      tree_ratio<-choose_grow_prune(proposal_table=proposal_table,proposal_matrix=proposal_matrix,prior_table=prior_tree[[i]],prior_matrix=prior_matrix[[i]],alpha=alpha,beta=beta)
    }
    
    if(node_grown==1){
      #get q(T'-> T)
      #get prob of going from prior tree to proposal tree which is
      #p(grow)*p(grow node chosen)*p(split var chosen)*p(split point chosen)
      
      #ratio<-choose_grow(ratio=ratio,tree_ratio=tree_ratio,prior_tree=list_of_tree_tables[[i]][[j-1]],proposal_tree=proposal_table)
      ratio<-choose_grow(ratio=ratio,tree_ratio=tree_ratio,prior_tree=prior_tree[[i]],proposal_tree=proposal_table)
      
    }
    if(node_pruned==1){
      #get probability of going from proposal tree to prior tree which is
      #p(prune)*p(prune node chosen)
      #ratio<-choose_prune(ratio=ratio,tree_ratio=tree_ratio,prior_tree=list_of_tree_tables[[i]][[j-1]],proposal_tree=proposal_table)
      ratio<-choose_prune(ratio=ratio,tree_ratio=tree_ratio,prior_tree=prior_tree[[i]],proposal_tree=proposal_table)
    }
    #now choose which tree to update
    #now update step with the most likely tree
    
    if(ratio<1){
      #if ratio<1 then set tree as proposal tree with probability ratio else keep the prior tree
      prob<-runif(1)
      #probs[i]<-prob
      if(prob<ratio){
        #if random number is less than ratio use proposal tree
        #post_prob[i]<-prop_lik
        # list_of_tree_tables[[i]][[j]]<-proposal_table
        #list_of_tree_matrix[[i]][[j]]<-proposal_matrix
        post_tree[[i]]<-proposal_table
        prior_tree[[i]]<-post_tree[[i]]
        post_matrix[[i]]<-proposal_matrix
        prior_matrix[[i]]<-post_matrix[[i]]
        #here I need to save the yhat and importance variables for each tree
        # terminal_nodes<-which(list_of_tree_tables[[i]][[j]][,5]==-1)
        terminal_nodes<-which(post_tree[[i]][,5]==-1)        
      }else{
        #else use prior
        
        # list_of_tree_tables[[i]][[j]]<- list_of_tree_tables[[i]][[j-1]]
        #  list_of_tree_matrix[[i]][[j]]<-list_of_tree_matrix[[i]][[j-1]]
        post_tree[[i]]<- prior_tree[[i]]
        
        if(nrow(prior_tree[[i]])==1){
          prior_matrix[[i]]<-as.matrix(prior_matrix[[i]])
        }
        post_matrix[[i]]<-prior_matrix[[i]]
      }
    }else{
      #if ratio>1 use proposal
      
      # list_of_tree_tables[[i]][[j]]<-proposal_table
      #list_of_tree_matrix[[i]][[j]]<-proposal_matrix
      #terminal_nodes<-which(list_of_tree_tables[[i]][[j]][,5]==-1)
      post_tree[[i]]<-proposal_table
      post_matrix[[i]]<-proposal_matrix
      prior_tree[[i]]<-post_tree[[i]]
      
      terminal_nodes<-which(post_tree[[i]][,5]==-1)
      if(nrow(prior_tree[[i]])==1){
        prior_matrix[[i]]<-as.matrix(prior_matrix[[i]])
      }
      prior_matrix[[i]]<-post_matrix[[i]]
      
    } 
    #if(class(list_of_tree_matrix[[i]][[j]])=="numeric"){
    if(class(post_matrix[[i]])=="numeric"){
      # list_of_tree_matrix[[i]][[j]]<-matrix(list_of_tree_matrix[[i]][[j]])
      post_matrix[[i]]<-matrix(post_matrix[[i]])
      colnames(proposal_table)<-c("left daughter","right daughter","split var","split point","status","mean","std dev")
    }
    
    #terminal_nodes<-which(list_of_tree_tables[[i]][[j]][,5]==-1)
    terminal_nodes<-which(post_tree[[i]][,5]==-1)
    
    # new_mean_var<-update_mean_var(tree_table=list_of_tree_tables[[i]][[j]],tree_matrix=list_of_tree_matrix[[i]][[j]],resids=resids,a=a)
    new_mean_var<-update_mean_var(tree_table=post_tree[[i]],tree_matrix=post_matrix[[i]],resids=resids,a=a,sigma=sigma,mu_mu=mu_mu)
    # set.seed(1)
    node_means<- rnorm(length(terminal_nodes),new_mean_var[[1]],sqrt(new_mean_var[[2]]))
    # node_means<- rnorm(length(terminal_nodes),new_mean_var2[[1]],sqrt(new_mean_var2[[2]]))
    
    # updated_preds<-update_predictions(tree_table=list_of_tree_tables[[i]][[j]],tree_matrix=list_of_tree_matrix[[i]][[j]],
    #                                  new_mean=node_means,new_var=new_mean_var[[2]],n=n)
    updated_preds<-update_predictions(tree_table=post_tree[[i]],tree_matrix=post_matrix[[i]],
                                      new_mean=node_means,new_var=new_mean_var[[2]],n=n)
    
    predy[,i]<-updated_preds[[2]]
    #list_of_tree_tables[[i]][[j]]<-updated_preds[[1]]
    post_tree[[i]]<-updated_preds[[1]]
    
    #now have chosen tree add the variable counts
    #vars chosen in posterior tree
    # var_imp[j,]<-get_imp_vars(tree_table=list_of_tree_tables[[i]][[j]][,3],num_col=ncol(data),current_vars=var_imp[j,])
    var_imp[j,]<-get_imp_vars(tree_table=post_tree[[i]][,3],num_col=ncol(data),current_vars=var_imp[j,])
    

    # REMOVE FROM INSIDE TREE LOOP --------------------------------------------

    #update tau for next iteration
    #set a and b (changing parameters from (nu.lambda) in inverse gamma to gamma)
    a1=nu/2
    b=2/(nu*lambda)
    
    #calculate residuals (response for tree)
    #sum_preds<-rep(0,n)
    
    #sum_preds<-apply(predy,1,sum)
    resids<-calculate_resids(predictions=predy,response=y_scaled)
    sum_preds<-calc_rowsums(predictions=predy,response=y_scaled)
    #resids_c<-calculate_resids(predictions=predy,response=y_scaled)
    #resids<-resids_c[[1]]
    #sum_preds<-resids_c[[2]]
    #update residuals for next iteration
    #resids<-y_scaled-sum_preds
    #s<-resids^2  
    #s<-sum(s)  
    #
    
    sigma<-update_sigma(a1=a1,b=b,resids=resids,n=n) 
    #sigma<-rgamma(1,shape=(a1+n/2),rate=((s/2)+(1/b)))  
    #sigma<-sqrt(1/sigma)  
    
    sigma_its[j]<-sigma 
    
    # REMOVE FROM INSIDE TREE LOOP --------------------------------------------
    
  }
  post_predictions[j,]<-apply(predy,1,sum)
  sum_preds<-apply(predy,1,sum)
  original_y<-get_original(low=min(y),high=max(y),sp_low=-0.5,sp_high=0.5,sum_preds=sum_preds)
  
  post_predictions_orig[j,]<-original_y
}
end=Sys.time()
start-end


#create a bart object and plot.bart function
par(mfrow=c(2,2))
plot(sigma_its,ylab="Sigma",col=c(rep(3,3000),rep(1,niter-3000)))
abline(h=sigma_true,col='red')

total_imp<-apply(var_imp,2,sum)
freq_imp<-total_imp/sum(total_imp)
quantile(freq_imp,0.9)
p=ncol(data)
#dev.off()
#par(mfrow=c(1,2))
plot(freq_imp,type="o",main="Variable Importance",xlim=c(1,p))
lines(c(1:p),rep(1/p,p))
##########################################################################
post_predictions_orig<-post_predictions_orig[-1,]

plquants=c(.05,.95)
cols =c('blue','black')
ql <- apply(post_predictions_orig[-c(1:3000),],2,quantile,probs=plquants[1])
qm <- apply(post_predictions_orig[-c(1:3000),],2,quantile,probs=.5)
qu <- apply(post_predictions_orig[-c(1:3000),],2,quantile,probs=plquants[2])
plot(y,qm,ylim=range(ql,qu),xlab='y',ylab=
       'Posterior Interval for E(Y|x)')
for (i in 1:length(qm))
  lines(rep(y[i],2),c(ql[i],qu[i]),col=cols[1])
abline(0,1,lty=2,col=cols[2])
