#include <Rcpp.h>
using namespace Rcpp;
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
  //  std::cout<<sv_col[obs_to_update[i]-1]<<"\n";
    
    if(sv_col[obs_to_update[i]]<=split_point){
    //  std::cout<<"xmat col2"<<sv_col[obs_to_update[i]-1]<<"\n";
     
      ld_obs.push_back(obs_to_update[i]);
      //std::cout<<"ld_obs"<<ld_obs[0]<<"\n";
    }
    if(sv_col[obs_to_update[i]]>=split_point){
    //  std::cout<<"xmat col2"<<sv_col[obs_to_update[i]-1]<<"\n";
     
      rd_obs.push_back(obs_to_update[i]);
      //std::cout<<"ld_obs"<<ld_obs[0]<<"\n";
    }
  }
  //ld_obs=xmat(obs_to_update[,1],)
  //ld_obs<-ld_obs[which(ld_obs[,split_var]<=split_point),]
  daughter_obs[0]=ld_obs;
  daughter_obs[1]=rd_obs;
  return(daughter_obs);
}

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


#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

NumericVector find_term_cols(IntegerMatrix tree_matrix,int terminal_node){
  NumericVector term_cols;
for(int i=0;i<tree_matrix.nrow();i++){
      for(int j=0;j<tree_matrix.ncol();j++){
        
        if(tree_matrix(i,j)==terminal_node){
          term_cols.push_back(j);  
        }
      }
    }
    return(term_cols);
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
IntegerMatrix set_daughter(int left_daughter,int right_daughter,IntegerVector ld_obs,IntegerVector rd_obs,IntegerMatrix tree_matrix,double term_cols){
  for(int i=0;i<ld_obs.size();i++){
    tree_matrix(ld_obs[i],term_cols+1)=left_daughter;    
  }
  for(int i=0;i<rd_obs.size();i++){
    tree_matrix(rd_obs[i],term_cols+1)=right_daughter;    
  }    
  return(tree_matrix);
}
//***************************************************************************************/
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

IntegerMatrix update_tree_matrix(NumericMatrix xmat,int swap_parent,NumericMatrix treedata,IntegerMatrix tree_matrix){
  int left_daughter;
  int right_daughter;
  left_daughter=treedata(swap_parent-1,0);
  right_daughter=treedata(swap_parent-1,1);
  int ind_nodes_to_update=1;
  IntegerVector nodes_to_update;
  nodes_to_update[0]=swap_parent;  
  int split_var;
  double split_point;
  IntegerVector ld_obs;
  IntegerVector rd_obs;
  int temp_ind_nodes_to_update=ind_nodes_to_update;
  IntegerVector temp_nodes_to_update=nodes_to_update;
  while(ind_nodes_to_update!=0){
  //go in here when there is still a node that was split so observations in children 
  //nodes need to be reassigned
  //for(int c=0;c<2;c++){
   // std::cout<<"ENTER LOOP"<<c<<"\n";
      //go in here when there is still a node that was split so observations in children 
      //nodes need to be reassigned
    for(int i=0;i<ind_nodes_to_update;i++){
    
      NumericVector obs_to_update=find_term_obs(tree_matrix,nodes_to_update[i]);
    
      NumericVector term_cols=unique(find_term_cols(tree_matrix,nodes_to_update[i]));
      //tree_matrix(obs_to_update,unique(find_term_cols(tree_matrix,nodes_to_update[i]))+1)=999;

    //  std::cout<<"obs_to_update"<<obs_to_update.size()<<"cols to update"<<term_cols[0]<<"\n";
    //  std::cout<<"left_daughter"<<left_daughter<<"right_daughter"<<right_daughter<<"\n";
      
      for(int j =0;j<obs_to_update.size();j++){
        tree_matrix(obs_to_update[j],term_cols[0]+1)=999;
      }
      
    //print(c("Check nodes to update",nodes_to_update))
      left_daughter=treedata(nodes_to_update[i]-1,0);    
      right_daughter=left_daughter+1;
      split_var=treedata(nodes_to_update[i]-1,2);
      split_point=treedata(nodes_to_update[i]-1,3);
      
      //this means that we haven't got to the end of the branch yet
        //follow it down and keep refining the subset until we get to a terminal node
      
      
      ld_obs=get_daughter_obs(xmat,obs_to_update,split_var,split_point)[0];
      rd_obs=get_daughter_obs(xmat,obs_to_update,split_var,split_point)[1];
     
     /*now have the new subsets which go to the left and right daughters of the change node
        need to update these in the tree_matrix
        we can only change internal nodes which means there will always be a child so will
        there will always be at least one column extra in the tree matrix for the children nodes*/
     std::cout<<"ld size "<<ld_obs[0]<<" "<<ld_obs[1]<<" "<<ld_obs[2]<<" "<<"rd_size "<<rd_obs[0]<<"term_cols[0]"<<term_cols[0]<<"\n";    
      tree_matrix=set_daughter(left_daughter,right_daughter,ld_obs,rd_obs,tree_matrix,term_cols[0]);  
      
      //if i>1 this is overwriting the nodes to update from previous internal nodes
        //at the same depth, try this
     // if(i==0){
      //  int temp_ind_nodes_to_update=ind_nodes_to_update;
      //  IntegerVector temp_nodes_to_update=nodes_to_update;
    //  }
     //now see if ld/rd have children nodes
        
    std::cout<<"ld "<<left_daughter<<"rd "<<right_daughter<<"sv "<<split_var<<"sp "<<split_point<<"\n";  
      //if i>1 this is overwriting the nodes to update from previous internal nodes
        //at the same depth, try this
    if(i==0){
      temp_ind_nodes_to_update=ind_nodes_to_update;
      temp_nodes_to_update=nodes_to_update;
      std::cout<<"These should be 1 "<<temp_ind_nodes_to_update<<" "<<temp_nodes_to_update[0]<<"\n";
    }
  std::cout<<"Left Daughter Depth is"<<term_cols[0]<<"total depth of tree"<< tree_matrix.ncol()<<"\n";
    //now see if ld/rd have children nodes
    if(treedata(left_daughter-1,4)==-1){
      std::cout<<"left daughter is terminal"<<"\n";
    /*when at the terminal node of the branch, change the treedata to update new node positions
    for the swapped tree.*/
      if(ld_obs.size()!=0){  
        
        if(term_cols[0]+2 < tree_matrix.ncol()){
          for(int j=(term_cols[0]+2);j<tree_matrix.ncol();j++){ 
            for(int k=0;k<ld_obs.size();k++){
              tree_matrix(ld_obs[k],j)=0;
            }
          }
        }
        
      }
    }else{
          temp_ind_nodes_to_update=temp_ind_nodes_to_update+1;
          temp_nodes_to_update.push_back(left_daughter);
          //left_daughter2<-treedata[left_daughter,1]
          //print("LD not terminal go again!")
          std::cout<<"ind nodes is"<<temp_ind_nodes_to_update<<"temp_nodes_to_update"<<temp_nodes_to_update[0]<<"\n";
    } 
    if(treedata(right_daughter-1,4)==-1){
          /*when at the terminal node of the branch, change the treedata to update the mean and 
          st dev of tree terminal nodes.*/
          if(rd_obs.size()!=0){
            /*treedata[right_daughter,6]<-mean(y[as.numeric(rownames(rd_obs))])
            treedata[right_daughter,7]<-sd(y[as.numeric(rownames(rd_obs))])
            
            if obs were reassigned and some which were previously non-terminal are now terminal
            then values for these obs in future columns have to be set to zero!*/
            
            if(term_cols[0]+2<tree_matrix.ncol()){
              for(int j=(term_cols[0]+2);j<tree_matrix.ncol();j++){ 
                for(int k=0;k<rd_obs.size();k++){
                  tree_matrix(rd_obs[k],j)=0;
                }
              }
              
            }
          }
        }else{
          temp_ind_nodes_to_update=temp_ind_nodes_to_update+1;
          temp_nodes_to_update.push_back(right_daughter); 
          std::cout<<"after rd ind nodes is"<<temp_ind_nodes_to_update<<"temp_nodes_to_update"<<temp_nodes_to_update[0]<<"size is"<<temp_nodes_to_update.size()<<"\n";
        }
  }
  // #need to check this when there is more than one child that needs to be updated
     // for(int k =0;k<ind_nodes_to_update;k++){
    //    temp_nodes_to_update.erase(k);        
    //  }
      IntegerVector temp;
      nodes_to_update=temp;
      nodes_to_update=temp_nodes_to_update;
      std::cout<<"nodes_to_update size"<<nodes_to_update.size()<<"\n";
      std::cout<<"nodes_to_update 1"<<nodes_to_update[0]<<"nodes_to_update 2"<<nodes_to_update[1]<<"\n";
      ind_nodes_to_update=temp_ind_nodes_to_update-ind_nodes_to_update;
      std::cout<<"ind_nodes size"<<ind_nodes_to_update<<"\n";
}    
             
  return(tree_matrix);
}