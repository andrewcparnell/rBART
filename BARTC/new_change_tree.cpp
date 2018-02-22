#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

NumericVector find_term_obs(IntegerMatrix tree_matrix_temp,int terminal_node){
  NumericVector term_obs;
  for(int i=0;i<tree_matrix_temp.nrow();i++){
      for(int j=0;j<tree_matrix_temp.ncol();j++){
        
        if(tree_matrix_temp(i,j)==terminal_node){
          term_obs.push_back(i);  
        }
      }
  }
  return(term_obs);
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

NumericVector find_term_cols(IntegerMatrix tree_matrix_temp,int terminal_node){
  NumericVector term_cols;
for(int i=0;i<tree_matrix_temp.nrow();i++){
      for(int j=0;j<tree_matrix_temp.ncol();j++){
        
        if(tree_matrix_temp(i,j)==terminal_node){
          term_cols.push_back(j);  
        }
      }
    }
    return(term_cols);
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
IntegerMatrix set_daughter(int left_daughter,int right_daughter,IntegerVector ld_obs,IntegerVector rd_obs,IntegerMatrix tree_matrix_temp,double term_cols){
  for(int i=0;i<ld_obs.size();i++){
    tree_matrix_temp(ld_obs[i],term_cols+1)=left_daughter;    
  }
  for(int i=0;i<rd_obs.size();i++){
    tree_matrix_temp(rd_obs[i],term_cols+1)=right_daughter;    
  }    
  return(tree_matrix_temp);
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

List get_daughter_obs(NumericMatrix xmat,NumericVector obs_to_update,int split_var,double split_point){
  IntegerVector ld_obs;
  IntegerVector rd_obs;
  NumericVector sv_col;
  List daughter_obs(2);
  //NumericVector sv_obs;
 
  sv_col=xmat(_,split_var-1);
 
  
  for(int i=0;i<obs_to_update.size();i++){
    
    if(sv_col[obs_to_update[i]]<=split_point){
 
     
      ld_obs.push_back(obs_to_update[i]);
 
    }
    if(sv_col[obs_to_update[i]]>=split_point){
 
     
      rd_obs.push_back(obs_to_update[i]);
 
    }
  }
 
  daughter_obs[0]=ld_obs;
  daughter_obs[1]=rd_obs;
  return(daughter_obs);
}

#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
List update_changetree_matrix(NumericMatrix xmat,NumericMatrix treedata_temp,IntegerMatrix tree_matrix_temp,double split_var,double split_point,int change_node){
  List ret(2);
  int left_daughter;
  int right_daughter;
  IntegerVector ld_obs;
  IntegerVector rd_obs;
  int ind_nodes_to_update=1;
  int temp_ind_nodes_to_update=ind_nodes_to_update;
  IntegerVector temp_nodes_to_update;
 //#assign new split var and point to the treedata
  treedata_temp(change_node-1,2)=split_var;
  treedata_temp(change_node-1,3)=split_point;
    
//#find left daughter and update tree_matrix and mean st dev if it is a terminal node
  left_daughter=treedata_temp(change_node-1,0);
  right_daughter=treedata_temp(change_node-1,1);
  
  
  IntegerVector nodes_to_update;
  nodes_to_update.push_back(change_node);
 while(ind_nodes_to_update!=0){ 
  for(int i=0;i<ind_nodes_to_update;i++){
    
    NumericVector obs_to_update=find_term_obs(tree_matrix_temp,nodes_to_update[i]);
    
    NumericVector term_cols=unique(find_term_cols(tree_matrix_temp,nodes_to_update[i]));
         
      
    for(int j =0;j<obs_to_update.size();j++){
      tree_matrix_temp(obs_to_update[j],term_cols[0]+1)=999;
    }
      
    
    left_daughter=treedata_temp(nodes_to_update[i]-1,0);    
    right_daughter=left_daughter+1;
    split_var=treedata_temp(nodes_to_update[i]-1,2);
    split_point=treedata_temp(nodes_to_update[i]-1,3);
    ld_obs=get_daughter_obs(xmat,obs_to_update,split_var,split_point)[0];
    rd_obs=get_daughter_obs(xmat,obs_to_update,split_var,split_point)[1];
     
    /*now have the new subsets which go to the left and right daughters of the change node
      need to update these in the tree_matrix
      we can only change internal nodes which means there will always be a child so will
      there will always be at least one column extra in the tree matrix for the children nodes*/
   
    tree_matrix_temp=set_daughter(left_daughter,right_daughter,ld_obs,rd_obs,tree_matrix_temp,term_cols[0]);  
        //now see if ld/rd have children nodes
     if(i==0){
      temp_ind_nodes_to_update=ind_nodes_to_update;
      //*************/
      //need to see what is inseide temp_ind_nodes here!!!
 
      IntegerVector temp_temp;
      temp_nodes_to_update=temp_temp;
      temp_nodes_to_update=nodes_to_update;
 
    }    
    if(treedata_temp(left_daughter-1,4)==-1){
 
    /*when at the terminal node of the branch, change the treedata to update new node positions
    for the swapped tree.*/
      if(ld_obs.size()!=0){  
        
        if(term_cols[0]+2 < tree_matrix_temp.ncol()){
          for(int j=(term_cols[0]+2);j<tree_matrix_temp.ncol();j++){ 
            for(int k=0;k<ld_obs.size();k++){
              tree_matrix_temp(ld_obs[k],j)=0;
            }
          }
        }
        
      }
    }else{
          temp_ind_nodes_to_update=temp_ind_nodes_to_update+1;
          temp_nodes_to_update.push_back(left_daughter);
          //left_daughter2<-treedata_temp[left_daughter,1]
          //print("LD not terminal go again!")
       
    }
    if(treedata_temp(right_daughter-1,4)==-1){
     
          /*when at the terminal node of the branch, change the treedata to update the mean and 
          st dev of tree terminal nodes.*/
          if(rd_obs.size()!=0){
            /*treedata[right_daughter,6]<-mean(y[as.numeric(rownames(rd_obs))])
            treedata[right_daughter,7]<-sd(y[as.numeric(rownames(rd_obs))])
            
            if obs were reassigned and some which were previously non-terminal are now terminal
            then values for these obs in future columns have to be set to zero!*/
            
            if(term_cols[0]+2<tree_matrix_temp.ncol()){
              for(int j=(term_cols[0]+2);j<tree_matrix_temp.ncol();j++){ 
                for(int k=0;k<rd_obs.size();k++){
                  tree_matrix_temp(rd_obs[k],j)=0;
                }
              }
              
            }
          }
        }else{
          temp_ind_nodes_to_update=temp_ind_nodes_to_update+1;
         
          temp_nodes_to_update.push_back(right_daughter); 
          
        }
  }  
  // #need to check this when there is more than one child that needs to be updated
 
     for(int k =0;k<ind_nodes_to_update;k++){
     
        temp_nodes_to_update.erase(k); 
        
      }
 
      IntegerVector temp;
      
      nodes_to_update=temp;
      nodes_to_update=temp_nodes_to_update;
     // temp_ind_nodes_to_update=0;
      temp_nodes_to_update=temp;
      
    
      ind_nodes_to_update=temp_ind_nodes_to_update-ind_nodes_to_update;
  
}    
  ret[0]=treedata_temp;
  ret[1]=tree_matrix_temp;
  return(ret);
}
//########################################################################################//
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]

List new_change_tree(NumericMatrix xmat,IntegerMatrix tree_matrix,NumericMatrix treedata,int split_var,double split_point,int change_node){
  NumericMatrix treedata_temp=clone(treedata);
  IntegerMatrix tree_matrix_temp=clone(tree_matrix);
  List ret(2);

  IntegerVector internal_nodes;

  if(treedata_temp.nrow()==1){
   
   ret[0]=treedata_temp;
   ret[1]=tree_matrix_temp;
  
  }else{
    ret=update_changetree_matrix(xmat,treedata_temp,tree_matrix_temp, split_var,split_point, change_node); 
  }
  return(ret);
}
