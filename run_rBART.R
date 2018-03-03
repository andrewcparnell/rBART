# Some code for Mahdi to run BART machine and return the posterior trees and precision values

# Clear the workspace and load in package
rm(list = ls())
library(bartMachine) # install.packages('bartMachine') if required
#source("https://bioconductor.org/biocLite.R")
library(treeio) # biocLite("treeio")
library(ggtree) # biocLite("ggtree")

# Set the working directory before running
# Load in the data - first column is target
data = read.table("friedman.txt")
X = data[,-1]
y = scale(data[,1])[,1]

# # Run it through bartMachine with a fixed seed and 1 tree
# set.seed(123)
# bart_machine = bartMachine(X, y, num_trees = 2)
# summary(bart_machine)
# 
# # Extract predictions and plot vs true values
# bart_machine_pred = bart_machine_get_posterior(bart_machine,
#                                                new_data = X)
# y_hat_bartm = bart_machine_pred$y_hat
# #y_hat_post = bart_machine_pred$y_hat_posterior_samples
# plot(y, y_hat_bartm) # Unsurprisingly pretty good
# 
# # Extract the sig squared posteriors
# sigsqs_bartm = get_sigsqs(bart_machine)
# sigma_bartm = sqrt(sigsqs_bartm)


source('rBART.R')
set.seed(123)
rBART_out = rBART(X, y, num_trees = 2, 
                  MCMC = list(iter = 10,
                              burn = 0,
                              thin = 1))

# Plot posterior sigma and compare with BARTMachine
# plot(rBART_out$sigma, ylim = range(c(rBART_out$sigma, sigma_bartm)))
# points(sigma_bartm, col = 'blue')

# Compare predictions from truth
y_hat_rBART = apply(rBART_out$y_hat, 2, 'mean')
plot(y, y_hat_rBART)
cor(y, y_hat_rBART)
abline(a = 0, b = 1)

# Compare predictions between rBART and BARTMachine
# plot(y_hat_bartm, y_hat_rBART)
# abline(a = 0, b = 1)
 
# Check predict function works
pred_y_rBART = predict_rBART(X, rBART_out, type = 'mean')
plot(pred_y_rBART, y_hat_rBART) # Should be identical
abline(a = 0, b = 1)


# Tree plotter ------------------------------------------------------------

source('rBART.R')
plot_tree(rBART_out, iter = 5, tree_num = 1, horiz = FALSE)
rBART_out$trees[[5]][[1]]$tree_matrix

