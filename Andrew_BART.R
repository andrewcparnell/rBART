# Create a version of BART that runs purely in R without any calls to other packages
# The idea is that this can be used for teaching purposes and so it clearly commented

# Boiler plate code -------------------------------------------------------


# Main function -----------------------------------------------------------

BART_Andrew = function(X, y, # X is the feature matrix, y is the target
                       control = list(num_trees = 2, # Number of trees
                                      node_min_size = 5), # Size of smallest nodes
                       priors = list(alpha = 0.95, # Prior control list
                                     beta = 2,
                                     mu_mu = 0,
                                     tau_mu = 2,
                                     nu = 3,
                                     lambda = 0.1), 
                       inits = list(tau = 1), # Initial values list
                       MCMC = list(iter = 1000, # Number of iterations
                                   burn = 200, # Size of burn in
                                   thin = 1) # Amount of thinning
                       ) { 
  
  # Extract control parameters
  num_trees = control$num_trees
  node_min_size = control$node_min_size
  
  # Extract hyper-parameters
  alpha = priors$alpha # Tree shape parameter 1
  beta = priors$beta # Tree shape parameter 2
  mu_mu = priors$mu_mu # Overall mean for terminal nodes
  tau_mu = priors$tau_mu # Precision for overall mean (sometimes called a?)
  nu = priors$nu # Parameter 1 for precision
  lambda = priors$lambda # Parameter 2 for precision
    
  # Extract initial values
  tau = inits$tau
  
  # Extract MCMC details
  iter = MCMC$iter # Number of iterations
  burn = MCMC$burn # Size of burn in
  thin = MCMC$thin # Amount of thinning
  
  # Storage containers
  store_size = (iter - burn)/thin
  tree_store = vector('list', store_size)
  sigma_store = rep(NA, store_size)
  y_hat_store = matrix(NA, ncol = length(y), nrow = store_size)

  # Scale the response target variable
  y_mean = mean(y)
  y_sd = sd(y)
  y_scale = (y - y_mean)/y_sd
  n = length(y_scale)
  
  # Create a list of trees for the initial stump
  curr_trees = create_stump(num_trees = num_trees, 
                            y = y_scale,
                            X = X)

  # Start the iterations loop
  for (i in 1:iter) {
    if(i%%10 ==0) cat('iteration',i,'\n')
    
    # If at the right place store everything
    if((i > burn) & ((i %% thin) == 0) ) {
      curr = (i - burn)/thin
      tree_store[[curr]] = curr_trees
      sigma_store[curr] = sigma
      y_hat_store[curr,] = predictions
    }
    
    # Start looping through trees
    for (j in 1:num_trees) {
      
      # Calculate partial residuals for current tree
      partial_trees = curr_trees 
      partial_trees[[j]] = NULL # Blank out that element of the list
      current_partial_residuals = y_scale - 
        get_predictions(partial_trees, X,
                        single_tree = num_trees == 2)
      
      # Propose a new tree via grow/change/prune/swap
      type = sample(c('grow', 'prune'), 1) #sample(c('grow', 'prune', 'change', 'swap'), 1)
      new_trees = curr_trees
      new_trees[[j]] = change_trees(y = y_scale,
                               X = X,
                               type = type, 
                               curr_tree = curr_trees[[j]],
                               node_min_size = node_min_size)

      # Check to see the tree is ok 
      if(any(new_trees[[j]]$tree_matrix[,'node_size'] == 0)) browser()
      
      # Calculate the complete conditional and acceptance probability
      new_log_lik = tree_full_conditional(new_trees[[j]], 
                                          current_partial_residuals,
                                          tau, 
                                          tau_mu)
      old_log_lik = tree_full_conditional(curr_trees[[j]], 
                                          current_partial_residuals,
                                          tau, 
                                          tau_mu)
      new_log_prior = get_tree_prior(new_trees[[j]], alpha, beta)
      old_log_prior = get_tree_prior(curr_trees[[j]], alpha, beta)

      # If accepting a new tree update all relevant parts
      accept_ratio = exp(new_log_lik + new_log_prior - old_log_lik - old_log_prior)
      
      # cat('iteration',i,'tree update',j,'\n')
      # cat('curr_tree\n')
      # print(curr_trees[[j]])
      # cat('new_tree\n')
      # print(new_trees[[j]])
      # cat('accept_ratio = ',accept_ratio,'\n')
      # browser()
      
      if(accept_ratio > runif(1)) {
        # Make changes if accept
        curr_trees = new_trees
        
        # Simulate mu for new trees
        curr_trees[[j]] = simulate_mu(new_trees[[j]], 
                                      current_partial_residuals, 
                                      tau, tau_mu)
        
      }
      
    } # End loop through trees
    
    # Calculate full set of predictions
    predictions = get_predictions(curr_trees, X)
    S = sum((y_scale - predictions)^2)
    
    # Update tau and sigma
    tau = update_tau(S, nu, lambda, 
                     n = length(y_scale))
    sigma = 1/sqrt(tau)
    #print(sigma)
    
  } # End iterations loop
  
  return(list(trees = tree_store,
         sigma = sigma_store*y_sd,
         y_hat = y_hat_store))
  
} # End main function


# Function to create stump ------------------------------------------------

create_stump = function(num_trees,
                        y,
                        X) {
  
  # Each tree has 6 columns and 2 elements
  # The 3 elements are the tree matrix, and the node indices
  # The tree matrix has columns:
  # Terminal (0 = no, 1 = yes)
  # Child left
  # Child right
  # Node parents
  # Split variable
  # Split value
  # mu values
  # Node size
  
  # Create holder for trees
  all_trees = vector('list', length = num_trees)
  # Loop through trees
  for (j in 1:num_trees) {
    # Set up each tree to have two elements in the list as described above
    all_trees[[j]] = vector('list', length = 2)
    # Give the elements names
    names(all_trees[[j]]) = c('tree_matrix', 
                              'node_indices')
    # Create the two elements: first is a matrix
    all_trees[[j]][[1]] = matrix(NA, ncol = 8, nrow = 1)
    
    # Second is the assignment to node indices
    all_trees[[j]][[2]] = rep(1, length(y))
    
    # Create column names
    colnames(all_trees[[j]][[1]]) = c('terminal', 
                                      'child_left', 
                                      'child_right',
                                      'parent',
                                      'split_variable',
                                      'split_value',
                                      'mu',
                                      'node_size')
    
    # Set values for stump 
    all_trees[[j]][[1]][1,] = c(1, 1, NA, NA, NA, NA, 0, length(y))
  } # End of loop through trees
  
  return(all_trees)
  
} # End of function

# Function to change trees ------------------------------------------------
# i.e. grow prune change swap 
  
change_trees = function(y, # Target variable
                        X, # Feature matrix
                        type = c('grow',  # Grow existing tree
                                 'prune', # Prune existing tree
                                 'change', # Change existing tree - change split variable and value for an internal node
                                 'swap'), # Swap existing tree - swap a parent/child combo where both are internal
                        curr_tree, # The current set of trees (not required if type is stump)
                        node_min_size) { # The minimum size of a node to grow
  
  # Each tree has 6 columns and 2 elements
  # The 3 elements are the tree matrix, and the node indices
  # The tree matrix has columns:
  # Terminal (0 = no, 1 = yes)
  # Child left
  # Child right
  # Node parents
  # Split variable
  # Split value
  # mu values
  # Node size

  # Create a copy of the curr_tree to change
  new_tree = curr_tree
  
  # The nod indices is just a vector displaying the terminal node number for each observation
  if(type == 'grow') {

    # Get the list of terminal nodes
    terminal_nodes = which(new_tree$tree_matrix[,'terminal'] == 1) # Create the list of terminal nodes
    
    # Find terminal node sizes
    terminal_node_size = new_tree$tree_matrix[terminal_nodes,'node_size']
    
    # Add two extra rows to the tree in question
    new_tree$tree_matrix = rbind(new_tree$tree_matrix,
                                  c(1, NA, NA, NA, NA, NA, NA, NA), # Make sure they're both terminal
                                  c(1, NA, NA, NA, NA, NA, NA, NA))
    
    # Choose a random terminal node to split 
    node_to_split = sample(terminal_nodes, 1, 
                           prob = as.integer(terminal_node_size > node_min_size)) # Choose which node to split, set prob to zero for any nodes that are too small
    
    # Choose a split variable uniformly from all columns
    split_variable = sample(1:ncol(X), 1)
    # Choose a split value from the range of the current node but stop it from choosing empty nodex
    low_bound = min(X[new_tree$node_indices == node_to_split,
                      split_variable]) + .Machine$double.eps
    high_bound = max(X[new_tree$node_indices == node_to_split,
                      split_variable]) - .Machine$double.eps
    split_value = runif(1, low_bound, high_bound)
    
    curr_parent = new_tree$tree_matrix[node_to_split, 'parent'] # Make sure to keep the current parent in there. Will be NA if at the root node
    new_tree$tree_matrix[node_to_split,1:6] = c(0, # Now not temrinal
                                                 nrow(new_tree$tree_matrix) - 1, # child_left is penultimate row
                                                 nrow(new_tree$tree_matrix),  # child_right is penultimate row
                                                 curr_parent,
                                                 split_variable,
                                                 split_value)
                                                             
    #  Fill in the parents of these two nodes
    new_tree$tree_matrix[nrow(new_tree$tree_matrix),'parent'] = node_to_split 
    new_tree$tree_matrix[nrow(new_tree$tree_matrix)-1,'parent'] = node_to_split 
    
    # Now call the fill function on this tree
    new_tree = fill_tree_details(new_tree, X)

  } # End of grow loop
  
  if(type == 'prune') {
    if(nrow(new_tree$tree_matrix) == 1) return(new_tree) # No point in pruning a stump!

    # Get the list of terminal nodes
    terminal_nodes = which(new_tree$tree_matrix[,'terminal'] == 1) # Create the list of terminal nodes
    
    # Pick a random termianl node to prune
    # ONLY PICK NODES WHERE BOTH LEFT AND RIGHT CHILD ARE TERMINAL
    bad_node_to_prune = TRUE # Assume a bad node pick
    while(bad_node_to_prune) {
      
      # Choose a random terminal node
      node_to_prune = sample(terminal_nodes, 1)
      
      # Find the parent of this terminal node
      parent_pick = new_tree$tree_matrix[node_to_prune, 'parent']
      
      # Get the two children of this parent
      child_left = new_tree$tree_matrix[parent_pick, 'child_left']
      child_right = new_tree$tree_matrix[parent_pick, 'child_right']
      
      # See whether either are terminal
      child_left_terminal = new_tree$tree_matrix[child_left, 'terminal']
      child_right_terminal = new_tree$tree_matrix[child_right, 'terminal']
      
      # If both are terminal then great
      if( (child_left_terminal == 1) & (child_right_terminal == 1) ) {
        bad_node_to_prune = FALSE # Have chosen a pair of terminal nodes so exist while loop
      }
        
    }# End of bad node to prune while loop
    
    # Delete these two rows from the tree matrix
    new_tree$tree_matrix = new_tree$tree_matrix[-c(child_left,child_right),,
                                                  drop = FALSE]
    # Make this node terminal again with no children or split values
    new_tree$tree_matrix[parent_pick,c('terminal',
                                        'child_left',
                                        'child_right',
                                        'split_variable',
                                        'split_value')] = c(1, NA, NA, NA, NA)
    
    # If we're back to a stump no need to call fill_tree_details
    if(nrow(new_tree$tree_matrix) == 1) {
      new_tree$node_indices = rep(1, length(y))
    } else {
      # If we've removed some nodes from the middle we need to re-number all the child_left and child_right values - the parent values will still be correct
      if(node_to_prune <= nrow(new_tree$tree_matrix)) { # Only need do this if we've removed some observations from the middle of the tree matrix
        # If you're pruning any nodes which affect parent indices further down the tree then make sure to shift the parent values
        bad_parents = which(new_tree$tree_matrix[,'parent']>=node_to_prune)
        # Shift them back because you have removed two rows
        new_tree$tree_matrix[bad_parents,'parent'] = new_tree$tree_matrix[bad_parents,'parent'] - 2

        for(j in node_to_prune:nrow(new_tree$tree_matrix)) {
          # Find the current parent
          curr_parent = new_tree$tree_matrix[j,'parent']
          # Find both the children of this node
          curr_children = which(new_tree$tree_matrix[,'parent'] == curr_parent)
          # Input these children back into the parent
          new_tree$tree_matrix[curr_parent,c('child_left','child_right')] = sort(curr_children)
        }
      }
      
      # Call the fill function on this tree
      new_tree = fill_tree_details(new_tree, X)
    }

  } # End of prune if statement
  
  # Return all_trees
  return(new_tree)
  
} # End of change_trees function

# Fill_tree_details -------------------------------------------------------

# The fill tree details function takes a tree matrix and returns the number of obs in each node in it and the indices of each observation in ecah terminal node
fill_tree_details = function(curr_tree, X) {
  
  # This code should essentially start from ignorance - no indices just a tree
  # Fill in the number of observations and the node indices
  
  # Collect right bits of tree
  tree_matrix = curr_tree$tree_matrix
  
  # Create a new tree matrix to overwrite
  new_tree_matrix = tree_matrix
  
  # Start with dummy node indices
  node_indices = rep(1, nrow(X))
  
  # For all but the top row, find the number of observations falling into each one
  for(i in 2:nrow(tree_matrix)) {
    # Get the parent
    curr_parent = tree_matrix[i,'parent']
    
    # Find the split variable and value of the parent
    split_var = tree_matrix[curr_parent,'split_variable']
    split_val = tree_matrix[curr_parent, 'split_value']
    
    # Find whether it's a left or right terminal node
    left_or_right = ifelse(tree_matrix[curr_parent,'child_left'] == i,
                           'left', 'right')
    if(left_or_right == 'left') {
      # If left use less than condition
      new_tree_matrix[i,'node_size'] = sum(X[node_indices == curr_parent,split_var] < split_val)
      node_indices[node_indices == curr_parent][X[node_indices == curr_parent,split_var] < split_val] = i
    } else {
      # If right use greater than condition
      new_tree_matrix[i,'node_size'] = sum(X[node_indices == curr_parent,split_var] >= split_val)
      node_indices[node_indices == curr_parent][X[node_indices == curr_parent,split_var] >= split_val] = i
    }
  } # ENd of loop through table
  return(list(tree_matrix = new_tree_matrix,
              node_indices = node_indices))  
  
} # End of function


# Get complete conditions -------------------------------------------------

tree_full_conditional = function(tree, R, tau, tau_mu) {
  # Function to compute log full conditional distirbution for an individual tree
  # R is a vector of partial residuals

  # Need to calculate log complete conditional, involves a sum over terminal nodes
  
  # First find which rows are terminal nodes
  which_terminal = which(tree$tree_matrix[,'terminal'] == 1)
  
  # Get node sizes for each terminal node
  nj = tree$tree_matrix[which_terminal,'node_size']
  
  # Get sum of residuals and sum of residuals squared within each terminal node
  sumRsq = aggregate(R, by = list(tree$node_indices), function(x) sum(x^2))[,2]
  sumR = aggregate(R, by = list(tree$node_indices), sum)[,2]
  
  # Now calculate the log posterior
  log_post = sum(0.5 * nj * log(tau) + 
    0.5 * log( tau_mu / (tau_mu + nj * tau)) -
    0.5 * tau * (sumRsq - tau * sumR^2 / (tau_mu + nj * tau) ) )
 
  return(log_post) 
}


# Get predictions ---------------------------------------------------------

# Gets the predicted values from a current set of trees
get_predictions = function(trees, X, single_tree = FALSE) {
  
  # Stop nesting problems in case of multiple trees
  if(is.null(names(trees)) & (length(trees) == 1)) trees = trees[[1]]
  
  # Normally trees will be a list of lists but just in case
  if(single_tree) {
    # Deal wiht just a single tree
    if(nrow(trees$tree_matrix) == 1) {
      predictions = rep(trees$tree_matrix[1, 'mu'], nrow(X))
    } else {
      # Loop through the node indices to get predictions
      predictions = rep(NA, nrow(X))
      unique_node_indices = unique(trees$node_indices)
      for(i in 1:length(unique_node_indices)) {
        predictions[trees$node_indices == unique_node_indices[i]] = 
          trees$tree_matrix[unique_node_indices[i], 'mu']
      }
    }
    # More here to deal with more complicated trees - i.e. multiple trees
  } else {
    # Do a recursive call to the function 
    partial_trees = trees 
    partial_trees[[1]] = NULL # Blank out that element of the list
    predictions = get_predictions(trees[[1]], X, single_tree = TRUE)  + 
      get_predictions(partial_trees, X, 
                      single_tree = length(partial_trees) == 1)
                      #single_tree = !is.null(names(partial_trees)))
    # The above only sets single_tree to if the names of the object is not null (i.e. is a list of lists)
  }
  
  return(predictions)
}


# Get tree priors ---------------------------------------------------------

get_tree_prior = function(tree, alpha, beta) {
  # Returns the tree log prior score

  # Need to work out the depth of the tree
  # First find the level of each node, then the depth is the maximum of the level
  level = rep(NA, nrow(tree$tree_matrix))
  level[1] = 1 # First row always level 1

  # Escpae quickly if tree is just a stump
  if(nrow(tree$tree_matrix) == 1) {
    return(log(alpha) - beta * log(2)) # Tree depth is 1 
  }
  
  for(i in 2:nrow(tree$tree_matrix)) {
    # Find the current parent
    curr_parent = tree$tree_matrix[i,'parent']
    # This child must have a level one greater than it's current parent
    level[i] = level[curr_parent] + 1
  }
  # Tree depth is max level
  tree_depth = max(level)
  
  # Calc and return log prior according to BART paper
  log_prior = log(alpha) - beta * log(1 + tree_depth)
  return(log_prior)
  
}

# Simulate_mu -------------------------------------------------------------

simulate_mu = function(tree, R, tau, tau_mu) {
  
  # Simulate mu values for a given tree
  
  # First find which rows are terminal nodes
  which_terminal = which(tree$tree_matrix[,'terminal'] == 1)
  
  # Get node sizes for each terminal node
  nj = tree$tree_matrix[which_terminal,'node_size']
  
  # Get sum of residuals in each terminal node
  sumR = aggregate(R, by = list(tree$node_indices), sum)[,2]
  
  # Now calculate mu values - NOTE THIS IS WRONG AND NEEDS ADJUSTING AND MATHS
  mu = rnorm(length(nj), 
             mean = tau * sumR / (nj * tau + tau_mu),
             sd = sqrt(1/(nj*tau + tau_mu)))
  
  # Wipe all the old mus out for other nodes
  tree$tree_matrix[,'mu'] = NA
  
  # Put in just the ones that are useful
  tree$tree_matrix[which_terminal,'mu'] = mu
  
  return(tree)
}

# Update tau --------------------------------------------------------------

update_tau = function(S, nu, lambda, n) {
  # THIS NEEDS TO BE FIXED BY DOING THE MATHS FIRST OF ALL
  tau = rgamma(1, (nu + n) / 2, (S + nu * lambda) / 2)
  
  return(tau)
}

