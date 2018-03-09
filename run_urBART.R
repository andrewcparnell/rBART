# This code produces a version of unsupervised Bayesian additive regression trees

# Clear the workspace and load in package
rm(list = ls())
library(treeio) # biocLite("treeio")
library(ggtree) # biocLite("ggtree")

# Source in code
source('urBART.R')

# Set the seed to sotp overfitting
set.seed(123)

# Load in the data - use the Old Faithful data
data(faithful)
y = scale(faithful)

urBART_out = urBART(X, y, num_trees = 2)

# Plot posterior sigma and compare with BARTMachine
plot(rBART_out$sigma, ylim = range(c(rBART_out$sigma, sigma_bartm)))
points(sigma_bartm, col = 'blue')
legend('topright', legend = c('rBART', 'BARTMachine'), pch = 1,
       col = c('black','blue'))

# Plot the log likelihood history
plot(rBART_out$log_lik)

# Compare predictions from truth
y_hat_rBART = apply(rBART_out$y_hat, 2, 'mean')
plot(y, y_hat_rBART)
cor(y, y_hat_rBART)
cor(y, y_hat_bartm) # Lower!
abline(a = 0, b = 1)

# Compare predictions between rBART and BARTMachine
plot(y_hat_bartm, y_hat_rBART)
abline(a = 0, b = 1)
 
# Check predict function works
pred_y_rBART = predict_rBART(X, rBART_out, type = 'mean')
plot(pred_y_rBART, y_hat_rBART) # Should be identical
abline(a = 0, b = 1)

# Tree plotter 
plot_tree(rBART_out, horiz = TRUE)

stop()

# WARNING: next part slow -------------------------------------------------

# 5-fold CV
rBART_out_CV = rBART_CV(X, y, num_trees = 10, folds = 5)
plot(y, rBART_out_CV$oob_predictions)
abline(a = 0, b = 1)
cor(y, rBART_out_CV$oob_predictions)

# 5-fold CV on BARTMachine
bartM_out_CV_pred = rep(NA, length(y))
fold_id = rBART_out_CV$fold_id
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
points(y, rBART_out_CV$oob_predictions, pch = 19)
abline(a = 0, b = 1)
cor(y,bartM_out_CV_pred)
