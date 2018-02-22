# Create a version of BART that runs purely in R without any calls to other packages
# The idea is that this can be used for teaching purposes and so it clearly commented

# Boiler plate code -------------------------------------------------------


# Main function -----------------------------------------------------------

BART_Andrew = function(X, y, # X is the feature matrix, y is the target
                       num_trees = 2, # Number of trees
                       priors = list(alpha = 0.95, # Prior control list
                                     beta = 2,
                                     mu_mu = 0,
                                     nu = 3,
                                     lambda = 0.1), 
                       inits = list(), # Initial values list
                       MCMC = list(iter = 100, # Number of iterations
                                   burn = 10, # Size of burn in
                                   thin = 1)) { # Amount of thinning
  
  # Extract hyper-parameters
  alpha = priors$alpha # Tree shape parameter 1
  beta = priors$beta # Tree shape parameter 2
  mu_mu = priors$mu_mu # Overall mean for terminal nodes
  nu = priors$nu # Parameter 1 for precision
  lambda = priors$lambda # Parameter 2 for precision
    
  # Extract initial values
  
  
  # Extract MCMC details
  iter = MCMC$iter # Number of iterations
  burn = MCMC$burn # Size of burn in
  thin = MCMC$thin # Amount of thinning
  
  # Scale the response target variable
  y_mean = mean(y)
  y_sd = sd(y)
  y_scale = (y - y_mean)/y_sd
  n = length(y_scale)
    
  # Create a list of trees for the initial stump
  curr_trees = create_trees(num_trees = num_trees, 
                            y = y_scale,
                            type = 'stump')
  browser()
  
  # Start the iterations loop
  for (i in 1:iter) {
    # Start looping through trees
    for (j in 1:num_trees) {
      
      # Calculate partial residuals for current tree
      current_partial_residuals = get_predictions(curr_trees[[-j]])
      
      # Propose a new tree via grow/change/prune/swap
      type = 'grow' #sample(c('grow', 'prune', 'change', 'swap'), 1)
      new_trees = create_trees(num_trees = num_trees,
                               y = y_scale,
                               type = type, 
                               which_tree = j, 
                               curr_trees = curr_trees)
      
      # Calculate the complete conditional and acceptance probability
      new_log_lik = get_log_lik(new_trees)
      old_log_lik = get_log_lik(curr_trees)
      new_log_prior = get_tree_prior(new_trees, j)
      old_log_prior = get_tree_prior(curr_trees, j)
      
      # If accepting a new tree update all relevant parts
      accept_ratio = exp(new_log_lik + new_log_prior - old_log_lik - old_log_prior)
      if(accept_ratio < runif(1)) {
        # Make changes if accept
      }
      
      
    } # End loop through trees
    
    # Calculate full set of predictions
    
    
    # Update sigma
    
    
    
  } # End iterations loop
  
  
  
  
} # End main function


# Function to create trees ------------------------------------------------

create_trees = function(num_trees, # Number of trees
                        y, # Target variable
                        type = c('stump', # Create initial stump trees
                                 'grow',  # Grow existing tree
                                 'prune', # Prune existing tree
                                 'change', # Change existing tree - change split variable and value for an internal node
                                 'swap'), # Swap existing tree - swap a parent/child combo where both are internal
                        which_tree, # Which tree to update (not required if type is stump)
                        curr_trees) { # The current set of trees (not required if type is stump)
  
  # Each tree has 10 columns and 2 elements
  # The 2 elements are the tree matrix and the node indices
  # The tree matrix has:
  # Node number
  # Terminal (0 = no, 1 = yes)
  # Child left
  # Child right
  # Node parents
  # Mean of target within node
  # mu values
  # Node size
  # Split variable
  # Split value
  
  # Create a stump
  if(type == 'stump') {
    # Create holder for trees
    all_trees = vector('list', length = num_trees)
    # Loop through trees
    for (j in 1:num_trees) {
      # Set up each tree to have two elements in the list as described above
      all_trees[[j]] = vector('list', length = 2)
      # Give the elements names
      names(all_trees[[j]]) = c('tree_matrix', 'node_indices')
      # Create the two elements: first is a matrix
      all_trees[[j]][[1]] = matrix(NA, ncol = 10, nrow = 1)
      
      # Second is the assignment to node indices
      all_trees[[j]][[2]] = rep(1, length(y))
      
      # Create column names
      colnames(all_trees[[j]][[1]]) = c('node', 
                                        'terminal', 
                                        'child_left', 
                                        'child_right',
                                        'parent',
                                        'y_mean',
                                        'mu',
                                        'size',
                                        'split_variable',
                                        'split_value')
      
      # Set values for stump 
      all_trees[[j]][[1]][1,] = c(1, 1, NA, NA, NA, 0, 0, length(y), NA, NA)
      
    } # End of loop through trees
    
  } # End of stump if statement

  if(type == 'grow') {
    browser()
    
    # Add two extra rows to the tree in question
    curr_trees[[which_tree]]$tree_matrix = rbind(curr_trees[[which_tree]]$tree_matrix, 0, 0)
    
    # Give these node numbers
    curr_trees[[which_tree]]$tree_matrix[,1] = 1:nrow(curr_trees[[which_tree]]$tree_matrix)
    
    # Choose a random terminal node to split 
    terminal_nodes = which(curr_trees[[which_tree]]$tree_matrix[,'terminal'] == 1) # Create the list of terminal ndoes
    node_to_split = sample(terminal_nodes, 1) # Choose which node to split
    curr_trees[[which_tree]]$tree_matrix[node_to_split,'terminal'] = 0 # Make that node non-terminal
    
    # 
    
    
  }
  
  # Return all_trees
  return(all_trees)
  
}


