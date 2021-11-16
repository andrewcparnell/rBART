# rBART2 uses varying sigma values for each terminal node
# See the associated maths.pdf document in the Github/maths folder for full details

# Clear the workspace and load in package
rm(list = ls())
# BiocManager::install(c("treeio", "ggtree"))
# library(treeio)
# library(ggtree)

# Source in the code
source('rBART2.R')

# Load in the data - first column is target
data = read.table("friedman.txt")
X = data[,-1]
y = scale(data[,1])[,1]
# Set the seed
set.seed(123)

# Simulate some data
# dat = sim_friedman(n = 200, p = 0)
# y = scale(dat$y)
# X = dat$X

rBART2_out = rBART2(X, y, num_trees = 10)
                    # MCMC = list(iter = 10,
                    #             burn = 0,
                    #             thin = 1))
stop()

# Compare predictions from truth
y_hat_rBART2 = apply(rBART2_out$y_hat, 2, 'mean')
plot(y, y_hat_rBART2)
abline(a = 0, b = 1)
cor(y, y_hat_rBART2)

# Check predict function works
pred_y_rBART2 = predict_rBART2(X, rBART2_out, type = 'mean')
plot(pred_y_rBART2, y_hat_rBART2) # Should be identical
abline(a = 0, b = 1)

# Tree plotter 
# plot_tree(rBART2_out, horiz = TRUE) # Doesn't work any more!

# Check uncertainties on the mean
pred_y_rBART2_mu = predict_rBART2(X, rBART2_out, type = 'all',
                                   par = 'mu')
# Now compute lower and upper 50% intervals
y_hat_rBART_mu_low_high = apply(rBART2_out$y_hat, 2, 'quantile',
                             c(0.25, 0.5, 0.75))
plot(y, y_hat_rBART_mu_low_high[2,])
sum <- 0
for(i in 1:ncol(y_hat_rBART_mu_low_high)) {
  lines(c(y[i], y[i]), 
        c(y_hat_rBART_mu_low_high[1,i], y_hat_rBART_mu_low_high[3,i]))
  if(y[i] < y_hat_rBART_mu_low_high[1,i] | y[i] > y_hat_rBART_mu_low_high[3,i]) {
    sum <- sum + 1
  }
}
title(1 - sum/ncol(y_hat_rBART_mu_low_high))
abline(a = 0, b = 1)

# Now do the same again but for tau. Should be larger than the previous one
pred_y_rBART2_tau = predict_rBART2(X, rBART2_out, type = 'all',
                                   par = 'tau')
# Now simulate from the model for these tau values for the above mu values
y_hat_full <- pred_y_rBART2_tau
for (i in 1:nrow(pred_y_rBART2_tau)) {
  for (j in 1:ncol(pred_y_rBART2_tau)) {
    y_hat_full[i,j] <- rnorm(1, pred_y_rBART2_mu[i,j], 
                             sd = 1/sqrt(pred_y_rBART2_tau[i,j]))
  }
}

y_hat_full_low_high = apply(y_hat_full, 2, 'quantile',
                             c(0.25, 0.5, 0.75))
plot(y, y_hat_full_low_high[2,]) # Should be identical
sum <- 0
for(i in 1:ncol(y_hat_full_low_high)) {
  lines(c(y[i], y[i]), 
        c(y_hat_full_low_high[1,i], y_hat_full_low_high[3,i]))
  if(y[i] < y_hat_full_low_high[1,i] | y[i] > y_hat_full_low_high[3,i]) {
    sum <- sum + 1
  }
}
abline(a = 0, b = 1)
title(1 - sum/ncol(y_hat_full_low_high))

# WARNING: next part slow -------------------------------------------------

# 5-fold CV
rBART2_out_CV = rBART2_CV(X, y, num_trees = 5, folds = 3)
plot(y, rBART2_out_CV$oob_predictions)
abline(a = 0, b = 1)
cor(y, rBART2_out_CV$oob_predictions)

