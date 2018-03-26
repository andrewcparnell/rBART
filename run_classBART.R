# Some code for Mahdi to run BART machine and return the posterior trees and precision values

# Clear the workspace and load in package
rm(list = ls())
library(bartMachine) # install.packages('bartMachine') if required
#source("https://bioconductor.org/biocLite.R")
library(treeio) # biocLite("treeio")
library(ggtree) # biocLite("ggtree")

# Source in the code
source('classBART.R')

# Load in the data - first column is target
data = read.table("friedman.txt")
X = data[,-1]
z = scale(data[,1])[,1]
y = rep(1, length(z))
y[z<0] = 0
# dat = sim_friedman(n = 200)
# y = as.vector(scale(dat$y))
# X = as.data.frame(dat$X)

# # Run it through bartMachine with a fixed seed and 1 tree
# set.seed(123)
# bart_machine = bartMachine(X, y, num_trees = 2)
# summary(bart_machine)
# # 
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

set.seed(123)
classBART_out = classBART(X, y, num_trees = 1,
                  MCMC = list(iter = 500,
                  burn = 0,
                  thin = 1))
y_hat_classBART = apply(classBART_out$y_hat, 2, 'mean')
plot(y, y_hat_classBART)
cor(y, y_hat_classBART)
stop()


# Plot the log likelihood history
plot(classBART_out$log_lik)

# Compare predictions from truth
y_hat_classBART = apply(classBART_out$y_hat, 2, 'mean')
plot(y, y_hat_classBART)
abline(a = 0, b = 1)
cor(y, y_hat_classBART)
cor(y, y_hat_bartm)

# Compare predictions between classBART and BARTMachine
plot(y_hat_bartm, y_hat_classBART)
abline(a = 0, b = 1)
 
# Check predict function works
pred_y_classBART = predict_classBART(X, classBART_out, type = 'mean')
plot(pred_y_classBART, y_hat_classBART) # Should be identical
abline(a = 0, b = 1)

# Tree plotter 
plot_tree(classBART_out, horiz = TRUE)

stop()

# WARNING: next part slow -------------------------------------------------

# 5-fold CV
classBART_out_CV = classBART_CV(X, y, num_trees = 10, folds = 5)
plot(y, classBART_out_CV$oob_predictions)
abline(a = 0, b = 1)
cor(y, classBART_out_CV$oob_predictions)

# 5-fold CV on BARTMachine
bartM_out_CV_pred = rep(NA, length(y))
fold_id = classBART_out_CV$fold_id
folds = 5
all_bartM_runs = vector('list', folds)
for (i in 1:folds) {
  cat('Running fold', i,'of',folds,'\n')
  X_in = X[fold_id !=i, , drop = FALSE]
  X_out = X[fold_id == i, , drop = FALSE]
  y_in = y[fold_id != i]
  y_out = y[fold_id ==i]
  
  all_bartM_runs[[i]] = bartMachine(X_in, y_in, 
                                    num_trees = 10)
  bart_machine_pred = bart_machine_get_posterior(bart_machine,
                                                 new_data = X_out)
  bartM_out_CV_pred[fold_id == i] = bart_machine_pred$y_hat
}

plot(y,bartM_out_CV_pred)
points(y, classBART_out_CV$oob_predictions, pch = 19)
abline(a = 0, b = 1)
cor(y,bartM_out_CV_pred)
