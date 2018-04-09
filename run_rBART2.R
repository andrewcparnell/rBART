# rBART2 uses varying sigma values for each terminal node
# See the associated maths.pdf document in the Github/maths folder for full details

# Clear the workspace and load in package
rm(list = ls())
#source("https://bioconductor.org/biocLite.R")
library(treeio) # biocLite("treeio")
library(ggtree) # biocLite("ggtree")

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
y_hat_rBART2 = apply(rBART2_out$y_hat, 2, 'mean')
plot(y, y_hat_rBART2)
abline(a = 0, b = 1)
cor(y, y_hat_rBART2)

#sqrt(mean(y - y_hat_rBART2)^2)

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
plot_tree(rBART2_out, horiz = TRUE)

stop()

# WARNING: next part slow -------------------------------------------------

# 5-fold CV
rBART2_out_CV = rBART2_CV(X, y, num_trees = 5, folds = 3)
plot(y, rBART2_out_CV$oob_predictions)
abline(a = 0, b = 1)
cor(y, rBART2_out_CV$oob_predictions)

