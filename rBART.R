# Create a version of BART that runs purely in R without any calls to other packages
# The idea is that this can be used for teaching purposes and so it clearly commented

# Main function -----------------------------------------------------------

rBART = function(X, y, # X is the feature matrix, y is the target
                 num_trees = 2, # Number of trees
                 control = list(node_min_size = 5), # Size of smallest nodes
                 priors = list(alpha = 0.95, # Prior control list
                               beta = 2,
                               tau_mu = 2,
                               nu = 3,
                               lambda = 0.1), 
                 inits = list(tau = 1), # Initial values list
                 MCMC = list(iter = 1250, # Number of iterations
                             burn = 250, # Size of burn in
                             thin = 1) # Amount of thinning
                 ) { 
  
  # Extract control parameters
  node_min_size = control$node_min_size
  
  # Extract hyper-parameters
  alpha = priors$alpha # Tree shape parameter 1
  beta = priors$beta # Tree shape parameter 2
  tau_mu = priors$tau_mu # Precision for overall mean (sometimes called a?)
  nu = priors$nu # Parameter 1 for precision
  lambda = priors$lambda # Parameter 2 for precision
    
  # Extract initial values
  tau = inits$tau
  sigma = 1/sqrt(tau)
  log_lik = 0

  # Extract MCMC details
  iter = MCMC$iter # Number of iterations
  burn = MCMC$burn # Size of burn in
  thin = MCMC$thin # Amount of thinning
  
  # Storage containers
  store_size = (iter - burn)/thin
  tree_store = vector('list', store_size)
  sigma_store = rep(NA, store_size)
  y_hat_store = matrix(NA, ncol = length(y), nrow = store_size)
  log_lik_store = rep(NA, store_size)
  
  # Scale the response target variable
  y_mean = mean(y)
  y_sd = sd(y)
  y_scale = (y - y_mean)/y_sd
  n = length(y_scale)
  
  # Create a list of trees for the initial stump
  curr_trees = create_stump(num_trees = num_trees, 
                            y = y_scale,
                            X = X)
  predictions = get_predictions(curr_trees, X, single_tree = num_trees == 1)
  
  # Set up a progress bar
  pb = utils::txtProgressBar(min = 1, max = iter,
                             style = 3, width = 60, 
                             title = 'Running rBART...')
  
  # Start the iterations loop
  for (i in 1:iter) {
    utils::setTxtProgressBar(pb, i)
    
    # If at the right place store everything
    if((i > burn) & ((i %% thin) == 0) ) {
      curr = (i - burn)/thin
      tree_store[[curr]] = curr_trees
      sigma_store[curr] = sigma
      y_hat_store[curr,] = predictions
      log_lik_store[curr] = log_lik
    }

    # Start looping through trees
    for (j in 1:num_trees) {
      
      # Calculate partial residuals for current tree
      if(num_trees > 1) {
        partial_trees = curr_trees 
        partial_trees[[j]] = NULL # Blank out that element of the list
        current_partial_residuals = y_scale - 
          get_predictions(partial_trees, X,
                          single_tree = num_trees == 2)
      } else {
        current_partial_residuals = y_scale
      }
      
      # Propose a new tree via grow/change/prune/swap
      type = sample(c('grow', 'prune', 'change', 'swap'), 1)
      if(i < max(floor(0.1*burn), 10)) type = 'grow' # Grow for the first few iterations 

      # Get a new tree!
      new_trees = curr_trees
      new_trees[[j]] = update_tree(y = y_scale,
                               X = X,
                               type = type, 
                               curr_tree = curr_trees[[j]],
                               node_min_size = node_min_size)

      # Calculate the complete conditional and acceptance probability
      l_new = tree_full_conditional(new_trees[[j]], 
                                          current_partial_residuals,
                                          tau, 
                                          tau_mu) + 
        get_tree_prior(new_trees[[j]], alpha, beta)
      
      l_old = tree_full_conditional(curr_trees[[j]], 
                                               current_partial_residuals,
                                               tau, 
                                               tau_mu) + 
        get_tree_prior(curr_trees[[j]], alpha, beta)
      
      # If accepting a new tree update all relevant parts
      a = exp(l_new - l_old)
      
      if(a > runif(1)) {
        # Make changes if accept
        curr_trees = new_trees
        
      } # End of accept if statement

      # Update mu whether tree accepted or not
      curr_trees[[j]] = simulate_mu(curr_trees[[j]], 
                                    current_partial_residuals, 
                                    tau, tau_mu)

    } # End loop through trees
    
    # Calculate full set of predictions
    predictions = get_predictions(curr_trees, X, single_tree = num_trees == 1)
    S = sum((y_scale - predictions)^2)
    
    # Update tau and sigma
    tau = update_tau(S, nu, lambda,
                     n = length(y_scale))
    sigma = 1/sqrt(tau)

    # Get the overall log likelihood
    log_lik = sum(dnorm(y_scale, mean = predictions, sd = sigma, log = TRUE))
    
  } # End iterations loop
  cat('\n') # Make sure progress bar ends on a new line
  
  return(list(trees = tree_store,
         sigma = sigma_store,
         y_hat = y_hat_store,
         log_lik = log_lik_store,
         y = y,
         X = X,
         iter = iter,
         burn = burn,
         thin = thin,
         store_size = store_size,
         num_trees = num_trees))
  
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

# Function to update trees ------------------------------------------------
# i.e. grow prune change swap 
  
update_tree = function(y, # Target variable
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

  # Call the appropriate function to get the new tree
  new_tree = switch(type,
                    grow = grow_tree(X, y, curr_tree, node_min_size),
                    prune = prune_tree(X, y, curr_tree),
                    change = change_tree(X, y, curr_tree),
                    swap = swap_tree(X, y, curr_tree))
  
  # Return the new tree
  return(new_tree)
  
} # End of update_tree function
  

# Grow_tree function ------------------------------------------------------

grow_tree = function(X, y, curr_tree, node_min_size) {
  
  # Set up holder for new tree
  new_tree = curr_tree
  
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
  # low_bound = min(X[new_tree$node_indices == node_to_split,
  #                   split_variable]) + .Machine$double.eps
  # high_bound = max(X[new_tree$node_indices == node_to_split,
  #                   split_variable]) - .Machine$double.eps
  # split_value = runif(1, low_bound, high_bound)
  
  # Alternatively follow BARTMachine and choose a split value using sample on the internal values of the available
  available_values = sort(unique(X[new_tree$node_indices == node_to_split,
                       split_variable]))
  split_value = sample(available_values[-c(1,length(available_values))], 1)
  
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

  # Return new_tree
  return(new_tree)
  
} # End of grow_tree function


# Prune_tree function -----------------------------------------------------

prune_tree = function(X, y, curr_tree) {
  
  # Create placeholder for new tree
  new_tree = curr_tree
  
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
      } # End for loop of correcting parents and children
    } # End if statement to fill in tree details
    
    # Call the fill function on this tree
    new_tree = fill_tree_details(new_tree, X)
    
  }
  
  # Return new_tree
  return(new_tree)
  
} # End of prune_tree function
  

# change_tree function ----------------------------------------------------

change_tree = function(X, y, curr_tree) {
  
  # Change a node means change out the split value and split variable of an internal node. Need to make sure that this does now produce a bad tree (i.e. zero terminal nodes)

  # If current tree is a stump nothing to change
  if(nrow(curr_tree$tree_matrix) == 1) return(curr_tree)
  
  # Create a holder for the new tree
  new_tree = curr_tree
  
  # Need to get the internal nodes
  internal_nodes = which(new_tree$tree_matrix[,'terminal'] == 0)
  terminal_nodes = which(new_tree$tree_matrix[,'terminal'] == 1)
  
  # Create a while loop to get good trees
  # Create a counter to stop after a certain number of bad trees
  max_bad_trees = 2
  count_bad_trees = 0
  bad_trees = TRUE
  while(bad_trees) {
    # Re-set the tree
    new_tree = curr_tree
    
    # choose an internal node to change
    node_to_change = sample(internal_nodes, 1)
    
    # Use the get_children function to get all the children of this node
    all_children = get_children(new_tree$tree_matrix, node_to_change)
    
    # Now find all the nodes which match these children
    use_node_indices = !is.na(match(new_tree$node_indices, all_children))
    
    # Create new split variable and value based on ignorance
    # then check this doesn't give a bad tree
    new_split_variable = sample(1:ncol(X), 1)
    available_values = sort(unique(X[use_node_indices,
                                     new_split_variable]))
    new_split_value = sample(available_values[-c(1,length(available_values))], 1)
    

    # Update the tree details
    new_tree$tree_matrix[node_to_change,
                         c('split_variable',
                           'split_value')] = c(new_split_variable, 
                                               new_split_value)
    
    # Update the tree node indices
    new_tree = fill_tree_details(new_tree, X)
    
    # Check for bad tree
    if(any(new_tree$tree_matrix[terminal_nodes, 'node_size'] == 0)) {
      count_bad_trees = count_bad_trees + 1
    } else {
      bad_trees = FALSE
    }
    if(count_bad_trees == max_bad_trees) return(curr_tree)
    
  } # end of while loop
  
  # Return new_tree
  return(new_tree)
  
} # End of change_tree function

# swap_tree function ------------------------------------------------------

swap_tree = function(X, y, curr_tree) {
  
  # Swap takes two neighbouring internal nodes and swaps around their split values and variables
  
  # If current tree is a stump nothing to change
  if(nrow(curr_tree$tree_matrix) == 1) return(curr_tree)
  
  # Create a holder for the new tree
  new_tree = curr_tree
  
  # Need to get the internal nodes
  internal_nodes = which(new_tree$tree_matrix[,'terminal'] == 0)
  terminal_nodes = which(new_tree$tree_matrix[,'terminal'] == 1)
  
  # If less than 3 internal nodes return curr_tree
  if(length(internal_nodes) < 3) return(curr_tree)
  
  # Find pairs of neighbouring internal nodes
  parent_of_internal = new_tree$tree_matrix[internal_nodes,'parent']
  pairs_of_internal = cbind(internal_nodes, parent_of_internal)[-1,]
  
  # Create a while loop to get good trees
  # Create a counter to stop after a certain number of bad trees
  max_bad_trees = 2
  count_bad_trees = 0
  bad_trees = TRUE
  while(bad_trees) {
    # Re-set the tree
    new_tree = curr_tree
    
    # Pick a random pair
    nodes_to_swap = sample(1:nrow(pairs_of_internal), 1)
    
    # Get the split variables and values for this pair
    swap_1_parts = new_tree$tree_matrix[pairs_of_internal[nodes_to_swap,1],
                                        c('split_variable', 'split_value')]
    swap_2_parts = new_tree$tree_matrix[pairs_of_internal[nodes_to_swap,2],
                                        c('split_variable', 'split_value')]
    
    # Update the tree details - swap them over
    new_tree$tree_matrix[pairs_of_internal[nodes_to_swap,1],
                         c('split_variable',
                           'split_value')] = swap_2_parts
    new_tree$tree_matrix[pairs_of_internal[nodes_to_swap,2],
                         c('split_variable',
                           'split_value')] = swap_1_parts
    
    # Update the tree node indices
    new_tree = fill_tree_details(new_tree, X)
    
    # Check for bad tree
    if(any(new_tree$tree_matrix[terminal_nodes, 'node_size'] == 0)) {
      count_bad_trees = count_bad_trees + 1
    } else {
      bad_trees = FALSE
    }
    if(count_bad_trees == max_bad_trees) return(curr_tree)
    
  } # end of while loop
  
  # Return new_tree
  return(new_tree)

} # End of swap_tree function

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
  } # End of loop through table
  
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
  sumRsq_j = aggregate(R, by = list(tree$node_indices), function(x) sum(x^2))[,2]
  S_j = aggregate(R, by = list(tree$node_indices), sum)[,2]

  # Now calculate the log posterior
  log_post = 0.5 * length(R) * log(tau) + sum(0.5 * log( tau_mu / (tau_mu + nj * tau)) -
    0.5 * tau * (sumRsq_j - tau * S_j^2 / (tau_mu + nj * tau) ) )
  return(log_post)
  # 
  # New Mahdi version - slower
  # P1 = 0.5 * length(R) * log(tau)
  # P2 = 0.5 * sum( log( tau_mu / (tau_mu + nj * tau)))
  # P3 = -0.5 * tau * sum( sumRsq_j )
  # P4 = 0.5 * (tau^2) * sum ( (S_j^2) / (tau_mu + nj * tau) )
  # 
  # return(P1 + P2 + P3 + P4)
}


# Get predictions ---------------------------------------------------------

# Gets the predicted values from a current set of trees
get_predictions = function(trees, X, single_tree = FALSE) {

  # Stop nesting problems in case of multiple trees
  if(is.null(names(trees)) & (length(trees) == 1)) trees = trees[[1]]
  
  # Normally trees will be a list of lists but just in case
  if(single_tree) {
    # Deal with just a single tree
    if(nrow(trees$tree_matrix) == 1) {
      predictions = rep(trees$tree_matrix[1, 'mu'], nrow(X))
    } else {
      # Loop through the node indices to get predictions
      predictions = rep(NA, nrow(X))
      unique_node_indices = unique(trees$node_indices)
      # Get the node indices for the current X matrix
      curr_X_node_indices = fill_tree_details(trees, X)$node_indices
      # Now loop through all node indices to fill in details
      for(i in 1:length(unique_node_indices)) {
        predictions[curr_X_node_indices == unique_node_indices[i]] = 
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
  level[1] = 0 # First row always level 0

  # Escpae quickly if tree is just a stump
  if(nrow(tree$tree_matrix) == 1) {
    return(log(1 - alpha)) # Tree depth is 0 
  }

  
  for(i in 2:nrow(tree$tree_matrix)) {
    # Find the current parent
    curr_parent = tree$tree_matrix[i,'parent']
    # This child must have a level one greater than it's current parent
    level[i] = level[curr_parent] + 1
  }
  
  # Only compute for the internal nodes
  internal_nodes = which(tree$tree_matrix[,'terminal'] == 0)
  log_prior = 0
  for(i in 1:length(internal_nodes)) {
    log_prior = log_prior + log(alpha) - beta * log(1 + level[internal_nodes[i]]) 
  } 
  # Now add on terminal nodes
  terminal_nodes = which(tree$tree_matrix[,'terminal'] == 1)
  for(i in 1:length(terminal_nodes)) {
    log_prior = log_prior + log(1 - alpha * ((1 + level[terminal_nodes[i]])^(-beta)))
  } 
  

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
  
  # Now calculate mu values
  mu = rnorm(length(nj) ,
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
  # Update from maths in Github folder
  tau = rgamma(1, shape = (nu + n) / 2, 
               rate = (S + nu * lambda) / 2)
  # Alternative
  #tau = rgamma(1, shape = (nu + n) / 2 - 1, scale = 2 / (S + nu * lambda))
  
  return(tau)
}


# get_children ------------------------------------------------------------

# A function which if, the current node is terminal, returns the node, or if not returns the children and calls the function again on the children
get_children = function(tree_mat, parent) {
  # Create a holder for the children
  all_children = NULL
  if(tree_mat[parent,'terminal'] == 1) {
    # If the node is terminal return the list so far
    return(c(all_children, parent))
  } else {
    # If not get the current children
    curr_child_left = tree_mat[parent, 'child_left']
    curr_child_right = tree_mat[parent, 'child_right']
    # Return the children and also the children of the children recursively
    return(c(all_children, 
             get_children(tree_mat,curr_child_left),
             get_children(tree_mat,curr_child_right)))
  }
}


# Predict function --------------------------------------------------------

predict_rBART = function(newX, rBART_posterior, 
                         type = c('all', 'median', 'mean')) {
  # Create predictions based on a new feature matrix
  # Note that there is minimal error checking in this - newX needs to be right!
  
  # Create holder for predicted values
  n_its = length(rBART_posterior$sigma)
  y_hat_mat = matrix(NA, nrow = n_its,
                     ncol = nrow(newX))
  
  # Now loop through iterations and get predictions
  for (i in 1:n_its) {
    # Get current set of trees
    curr_trees = rBART_posterior$trees[[i]]
    # Use get_predictions function to get predictions
    y_hat_mat[i,] = get_predictions(curr_trees, 
                                    newX, 
                                    single_tree = length(curr_trees) == 1)
  }

  # Sort out what to return
  out = switch(type,
               all = y_hat_mat,
               mean = apply(y_hat_mat,2,'mean'),
               median = apply(y_hat_mat,2,'median'))
  
  return(out)
  
} # end of predict function


# Plot_tree function ------------------------------------------------------

plot_tree = function(rBART_posterior, 
                     iter = sample(1:rBART_posterior$store_size, 1), 
                     tree_num = sample(1:rBART_posterior$num_trees, 1),
                     horiz = TRUE) {
  # Use the ggtree packages - needs to be installed from bioconductor
  require(ggtree)
  
  # cat out the iteration and tree 
  cat('iteraton:',iter,' tree number', tree_num,'\n')
  
  # Get the current tree we want to plot
  curr_tree = rBART_posterior$trees[[iter]][[tree_num]]$tree_matrix
  
  # Get rid of stumps
  if(nrow(curr_tree) ==1) stop("Tree is just a stump; cannot be plotted.")
  
  # Now need to get it into phylo format for plotting in ggtree
  
  # In a phylo tree the terminal nodes are labelled 1, 2, 3, etc. Don't know why!
  terminal_nodes = which(curr_tree[,'terminal'] == 1)
  internal_nodes = which(curr_tree[,'terminal'] == 0)
  
  # Create a vector which includes the new node numbers next to the old ones. The index is the old row number
  node_number_swap = cbind(c(terminal_nodes,internal_nodes),
                           1:nrow(curr_tree))
  node_number_vec = node_number_swap[order(node_number_swap[,1]),2]
  
  # This means I need to re-label the entire tree!
  new_tree = curr_tree[c(terminal_nodes, internal_nodes),]
  new_tree[,'child_left'] = node_number_vec[new_tree[,'child_left']]
  new_tree[,'child_right'] = node_number_vec[new_tree[,'child_right']]
  new_tree[,'parent'] = node_number_vec[new_tree[,'parent']]
  
  # Create the edge labels for the new tree
  split_variables = new_tree[,'split_variable']
  split_values = round(new_tree[,'split_value'],2)
  left_edge_labels = paste0('x',split_variables,'<',split_values)
  left_edge_labels[is.na(split_variables)] = NA
  right_edge_labels = paste0('x',split_variables,'>=',split_values)
  right_edge_labels[is.na(split_variables)] = NA
  
  # Get the new terminal nodes and internal nodes
  new_terminal_nodes = which(new_tree[,'terminal'] == 1)
  new_internal_nodes = which(new_tree[,'terminal'] == 0)
  
  # Now create the edge matrix needs to link each edge
  # Create holder for edge matrix
  edge = data.frame(parent = rep(NA,nrow(new_tree)*2), 
                    child = rep(NA,nrow(new_tree)*2), 
                    label = rep('x', nrow(new_tree*2)),
                    stringsAsFactors = FALSE)
  for(i in 1:nrow(new_tree)) {
      edge[2*i-1,1:2] = as.integer(c(new_tree[i,'child_left'],i))
      edge[2*i,1:2] = as.integer(c(new_tree[i,'child_right'],i))
      edge[2*i-1,3] = left_edge_labels[i]
      edge[2*i,3] = right_edge_labels[i]
  }
  # Trim out the NAs
  bad_rows = which(is.na(edge[,1]))
  tree_obj = edge[-bad_rows,c(2,1,3)]
  
  # Re-order it so that it follows the bizarre ordering of ape/ggtree
  tree_obj2 = tree_obj[1,1:3,drop = FALSE]
  tree_remaining = tree_obj[-1,,drop = FALSE]
  still_going = TRUE
  while(still_going) {
    curr_terminal = tree_obj2[nrow(tree_obj2),2]
    matches = which(tree_remaining[,1] == curr_terminal)
    if(length(matches) > 0) {
      tree_obj2 = rbind(tree_obj2, tree_remaining[matches,])  
      tree_remaining = tree_remaining[-matches,,drop=FALSE]
    } else {
      tree_obj2 = rbind(tree_obj2, tree_remaining[1,])  
      tree_remaining = tree_remaining[-1,,drop=FALSE]
    }
    if(nrow(tree_obj2) == nrow(tree_obj)) still_going = FALSE
  }
  
  # Get the number of internal nodes + stump node
  Nnode = length(new_internal_nodes)
  
  # Get the tip labels and the internal node labels
  tip_labels = round(new_tree[new_terminal_nodes, 'mu'], 2)
  
  # Put together 
  tree = list(edge = as.matrix(tree_obj2[,1:2]),
              tip.label = as.character(tip_labels),
              Nnode = Nnode)
  class(tree) = "phylo"
  attr(tree, 'order') = 'cladewise'
  
  # Now plot
  p = ggtree(tree) + 
    geom_tiplab(align = TRUE) + labs(title=paste0("Iteration: ",iter,', tree:',tree_num))
  if(!horiz) p = p + coord_flip() + 
    scale_x_reverse()
    
  # Add in edge labels
  edge = tree_obj2
  colnames(edge)=c("parent", "node", "edge_lab")
  
  # Finally combine and print
  p2 = p %<+% edge + geom_label(aes(x=branch, label=edge_lab))
  suppressWarnings(print(p2))

}


# rBART_CV ----------------------------------------------------------------

rBART_CV = function(X, y, folds = 5, num_trees = 2, ...) {
  # Run rBART with k-fold cross validation
  fold_id = sample(1:folds, size = length(y), replace=TRUE)
  
  # Create holder for predictions
  oob_predictions = rep(NA, length(y))
  # Holder for trees
  all_rBART_runs = vector('list', folds)
  
  # Now loop through and create new predictions
  for (i in 1:folds) {
    cat('Running fold', i,'of',folds,'\n')
    X_in = X[fold_id !=i, , drop = FALSE]
    X_out = X[fold_id == i, , drop = FALSE]
    y_in = y[fold_id != i]
    y_out = y[fold_id ==i]
    
    all_rBART_runs[[i]] = rBART(X_in, y_in, 
                                num_trees = num_trees, ...)
    oob_predictions[fold_id == i] = predict_rBART(newX = X_out, 
                                                  all_rBART_runs[[i]],
                                                  type = 'mean')
  }

  # Return output
  return(list(all_rBART_runs = all_rBART_runs,
              oob_predictions = oob_predictions,
              fold_id = fold_id))
}

# Simulate friedman -------------------------------------------------------

sim_friedman = function(n, p = 0, d = 1, scale_par = 5, scale_err = 0.5) {
  # Simulate some data using a multivariate version of Friedman
  # y = 10sin(πx1x2)+20(x3−0.5)2+10x4+5x5+ε
  X = matrix(NA, nrow = n, ncol = 5 + p)
  for(i in 1:ncol(X)) X[,i] = rnorm(n, 0, 1)
  pars = matrix(rnorm(5 * d, sd = scale_par), ncol = d, nrow = 5) # 5 parameters on d dimensions
  y = mean = matrix(NA, ncol = d, nrow = n)
  Sigma = rWishart(1, d, scale_err*diag(d))[,,1]
  if(d > 1) {
    err = rmvnorm(n, sigma = Sigma)
  } else {
    err = matrix(rnorm(n, sd = sqrt(Sigma)), ncol = 1)
  }
  for(j in 1:d) {
    mean[,j] = pars[1,j]*sin(X[,1]*X[,2]) + pars[2,j] * (X[,3]-0.5)^2 + pars[3,j] * X[,4] + pars[5,j] * X[,5]
    y[,j] = mean[,j] + err[,j]
  }
  return(list(y = y, X = X, Sigma = Sigma, mean = mean))
}

