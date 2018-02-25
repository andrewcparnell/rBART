# Some code for Mahdi to run BART machine and return the posterior trees and precision values

# Clear the workspace and load in package
rm(list = ls())
library(bartMachine) # install.packages('bartMachine') if required

# Set the working directory before running
# Load in the data - first column is target
data = read.table("friedman.txt")
X = data[,-1]
y = scale(data[,1])[,1]

# # Run it through bartMachine with a fixed seed and 1 tree
set.seed(123)
bart_machine = bartMachine(X, y, num_trees = 2)
summary(bart_machine)

# Extract predictions and plot vs true values
bart_machine_pred = bart_machine_get_posterior(bart_machine,
                                               new_data = X)
y_hat = bart_machine_pred$y_hat
#y_hat_post = bart_machine_pred$y_hat_posterior_samples
plot(y, y_hat) # Unsurprisingly pretty good

# Extract the sig squared posteriors
sigsqs = get_sigsqs(bart_machine)
sigma = sqrt(sigsqs)
plot(sigma)
# 
# output = cbind(tau, t(y_hat_post))
# colnames(output) = c('tau', paste0('y_hat',1:nrow(data)))
# 
# write.table(output, file = "bartmachine_output_t50.txt", quote = FALSE, row.names = FALSE)
# 
# # Run through Belinda's version
# # Rcpp::sourceCpp('BART_main_c++.cpp')

source('Andrew_BART.R')
set.seed(124)
BART_out = BART_Andrew(X, y, control = list(num_trees = 20, node_min_size = 5))
plot(BART_out$sigma)
y_hat_AP_mean = apply(BART_out$y_hat, 2, 'mean')
plot(y, y_hat_AP_mean)
abline(a = 0, b = 1)

plot(y_hat, y_hat_AP_mean)
abline(a = 0, b = 1)
 