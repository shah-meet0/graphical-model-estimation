# Ideas adapted from Meinshausen and Buhlmann (2006),
# Friedman, Hastie and Tibshirani (2007)
# Author: Meet Shah
# See Report for full list of references.


# install.packages('MASS')
# install.packages('glmnet')
# install.packages('glasso')
# install.packages('MESS)
# install.packages('CVglasso)

# INITIALIZATION (from here to end initialization there are only function definitions)

library(Matrix) 
library(MASS) # for mvrnorm
library(glmnet) # for implementing node-wise lasso
library(glasso) # for implementing glasso
library(MESS) # for Area under curve
library(CVglasso) # for cross-validating glasso Galloway (2018)

# MAKING SAMPLE + INFERENCE

# Whenever encountered:
# n - sample size
# p - number of predictors
# s_t - sparsity threshold, probability of conditional independence between two nodes.

# Constructs Theta  + Delta * Identity as called for in project.
covariance_matrix_maker <- function(predictors, sparsity_threshold, delta, uniform = FALSE){
  unif_matrix = matrix(runif(predictors^2), ncol=predictors)
  unif_matrix[lower.tri(unif_matrix)] = t(unif_matrix)[lower.tri(unif_matrix)]
  beta_matrix = matrix(rep(0, predictors^2), ncol = predictors)
  if(uniform){
    beta_matrix[unif_matrix >= sparsity_threshold] = runif(1,-1,1)
  }
  else{
    beta_matrix[unif_matrix >= sparsity_threshold] = 0.5
  }
  diag(beta_matrix) = delta
  return(beta_matrix)
}

# Standardizes to have unit diagonal. General purpose function,
# could have simply divided by delta.

covariance_matrix_standardizer <- function(covariance_matrix){
  D = diag(sqrt(diag(covariance_matrix)))
  Dinv = diag(1/diag(D))
  correlation_matrix = Dinv %*% covariance_matrix %*% Dinv
  return(correlation_matrix)
}


# Uses the previous functions to give a precision matrix with unit diagonals

get_correlation_matrix <- function(predictors, sparsity_threshold, delta, uniform = FALSE){
  A = covariance_matrix_maker(predictors, sparsity_threshold, delta, uniform)
  corr_matrix = covariance_matrix_standardizer(A)
  if(sum(eigen(corr_matrix)$values <= 0.00000001) > 0){
    return(get_correlation_matrix(predictors, sparsity_threshold, delta+0.1))
  }
  return(corr_matrix)
}


generate_multivariate_norm_sample <- function(num_observation, covar_matrix){
  cov_solved = chol2inv(chol(covar_matrix)) # Uses the fact that the precision matrix is symmetric positive definite
  p = length(covar_matrix[1,])
  mu = rep(0, p)
  sample = mvrnorm(n=num_observation, mu = mu, Sigma = cov_solved)
  return(sample)
}

lasso_sample_maker <- function(num_observation, predictors, sparsity_threshold,
                               delta, uniform = FALSE){
  theta = covariance_matrix_maker(predictors, sparsity_threshold, delta, uniform)
  std_theta = covariance_matrix_standardizer(theta)
  eigenvalues = eigen(std_theta)$values
  if(sum(eigenvalues <= 0.00000001) > 0){
    return(lasso_sample_maker(num_observation, predictors, sparsity_threshold, 
                              delta+0.1, uniform))
  }
  sample = generate_multivariate_norm_sample(num_observation, std_theta)
  return(list(sparsity_pattern = std_theta, sample = sample))
}

# Looks at the true edge set, and estimated edge set (Symmetric Matrices)
# And returns True Positive Rate, False Positive Rate, True Negative Rate
# False Negative rate, Proportion of correct edge predictions, 
# Misclassification Error rate
# Call with diagonal = True if the estimated edge set has diagonals set to True

summary_stat <- function(true_edge_set, edge_set, diagonal = FALSE){
  if(diagonal){
    diag(edge_set) = FALSE
  }
  n_params = sqrt(length(true_edge_set))
  # Divided by two to avoid double counting, as the edge sets formed 
  # count each edge twice.
  # If using edge sets, remember to take into account diagonals 
  # and double counting.
  true_positive = sum(true_edge_set == TRUE & edge_set == TRUE)/2
  false_positive = sum(true_edge_set == FALSE & edge_set == TRUE)/2
  true_negative = sum(true_edge_set == FALSE & edge_set == FALSE)/2
  false_negative = (sum(true_edge_set == TRUE & edge_set == FALSE) - n_params)/2
  
  tot_negative = sum(true_edge_set == FALSE)/2
  tot_positive = (sum(true_edge_set) - n_params)/2
  total_predictions = tot_negative + tot_positive
  total_correct = sum(true_edge_set == edge_set)/2/total_predictions 
  total_wrong = (sum(true_edge_set != edge_set) - n_params)/2/total_predictions
  
  true_positive_rate = true_positive/tot_positive
  false_positive_rate = false_positive/tot_negative
  true_negative_rate = true_negative/tot_negative
  false_negative_rate = false_negative/tot_positive
  
  summary_results = c(true_positive_rate, false_positive_rate, true_negative_rate, 
                      false_negative_rate, total_correct, total_wrong)
  summary_results = matrix(summary_results, nrow = 1)
  
  colnames(summary_results) <- c('TPR', 'FPR', 'TNR', 'FNR', 'TC %', 'Error %')
  # Note that actually TC % is the number reported times 100, 
  # same with Error %.
  # They are not reported in % to be on same scale as rest of the values.
  return(summary_results)
}

# NODEWISE LASSO

# Constructs MBpen lambda 
make_nodewise_lambda <- function(n, p, variable_var, alpha = 0.05){
  prob = alpha/(2 * p^2)
  cdf = pnorm(prob, lower.tail = FALSE)
  lambda  = 2 * sqrt(variable_var/n) * cdf
  return(lambda)
}

# Boolean and controls the way in which the edge set is formed.
nodewise_edgeset <- function(beta_matrix, and= TRUE, p = length(beta_matrix[,1])) {
  relation = beta_matrix != 0
  edge_set = matrix(0,ncol = p, nrow = p)
  if(and){
    edge_set = (relation == 1) & (t(relation) == 1)
  }
  else {
    edge_set = (relation == 1) | (t(relation) == 1)
  }
  
  return(edge_set)
}

# Primary function for performing node-wise lasso.
# method - MBpen uses a single value for lambda, Cv does k fold cross validation,
# lambdalist takes as input a list of lambda and returns an edge set for each of them
# alpha - controls MBpen value
# lambdalist - list of lambda provided if method = 'lambdalist'
# k - controls number of folds

# Returns edge set(s), fits of the function and lambdas used.

ndw_lasso <- function(data, alpha = 0.05, k = 5, lambdalist = c(0.01, 0.1, 1), and = FALSE, 
                      method = c("mbpen", "cv", "lambdalist")) {
  n = length(data[,1])
  p = length(data[1,])
  betas = matrix(0, nrow= p, ncol = p)
  lambdas = rep(0, p)
  fits  = list()
  method = method[1]
  if(method == "mbpen"){
    data = scale(data) # scale first to avoid ambiguity when choosing mbpen variance
    for(i in 1:p){
      lambdas[i] = make_nodewise_lambda(n, p, 1, alpha)
      glmfit = glmnet(x = data[, -i], y = data[,i], family = "gaussian", lambda = lambdas[i], intercept= FALSE, standardize = FALSE)
      betas[-i, i] = matrix(glmfit$beta)
      fits[[i]] = glmfit
    }
    edge_set = nodewise_edgeset(betas, and, p)
    return(list(edge_set = edge_set, fits = fits, lambdas = lambdas))
  }
  
  else if(method == "cv") {
    cv_fits = list()
    for(i in 1:p){
      fit_cv = cv.glmnet(x = data[,-i], y= data[,i], nfolds = k, family  = "gaussian", intercept = FALSE)
      lambdas[i] = fit_cv$lambda.min
      betas[-i, i] = matrix(coef(fit_cv))[-1]
      fits[[i]] = fit_cv
    }
    edge_set = nodewise_edgeset(betas, and, p)
    return(list(edge_set = edge_set, fits = fits, lambdas = lambdas))
  }
  
  else {
    lambda_fits = list()
    edge_sets = list()
    n_lambda = length(lambdalist)
    edge_seq = seq(from = 0, to = (n_lambda - 1) * p, by = p)
    betas = matrix(0, nrow = p * n_lambda, ncol= p)
    
    for(i in 1:p){
      i_seq = edge_seq + i
      glmfit = glmnet(x = data[, -i], y = data[,i], family = "gaussian", lambda = lambdalist, intercept = FALSE)
      lambda_fits[[i]] = glmfit
      betas[-i_seq, i] = matrix(glmfit$beta, ncol=1) 
    }
    
    for(j in 1:n_lambda) {
      if(j == n_lambda){
        row_seq = (edge_seq[j] + 1) : length(betas[,1])
      }
      else {
        row_seq = (edge_seq[j] + 1) : edge_seq[j+1]
      }
      beta_matrix_to_use = betas[row_seq, ]
      edge_sets[[j]] = nodewise_edgeset(beta_matrix_to_use, and, p)
    }
    return(list(edge_set = edge_sets, fits = lambda_fits, lambdas=lambdalist))
  }  
}

# Repeats the nodewise lasso multiple times given n,p, s_t and delta
# ... <- additional arguments to pass to ndwlasso function
ndw_lasso_repeater <- function(n, p , s_t, delta, times_repeated= 50, ...){
  results = matrix(rep(0, times_repeated * 6), nrow = times_repeated)
  sparsity_pattern = get_correlation_matrix(p, s_t, delta)
  true_edge_set = sparsity_pattern != 0
  progress_bar = txtProgressBar(min = 0, max = times_repeated, style =3)
  for (repitition in 1:times_repeated){
    generated_sample = generate_multivariate_norm_sample(n, sparsity_pattern)
    ndwresult = ndw_lasso(generated_sample, ...)
    summary = summary_stat(true_edge_set, ndwresult$edge_set)
    results[repitition, ] = summary
    setTxtProgressBar(progress_bar, repitition)
  }
  close(progress_bar)
  colnames(results) = c('TPR', 'FPR', 'TNR', 'FNR', 'TC %', 'Error %')
  return(results)
}

ndw_lasso_over_n <- function(list_n, p, s_t, delta, ...){
  sparsity_pattern = get_correlation_matrix(p, s_t, delta)
  true_edge_set = sparsity_pattern != 0
  n_iterations = length(list_n)
  results = matrix(rep(0, n_iterations*6), nrow = n_iterations)
  
  
  for(i in 1:n_iterations){
    sample = generate_multivariate_norm_sample(list_n[i], sparsity_pattern)
    edge_set = ndw_lasso(sample, ...)$edge_set
    results[i, ]= summary_stat(true_edge_set, edge_set)  
  }
  
  colnames(results) = c('TPR', 'FPR', 'TNR', 'FNR', 'TC %', 'Error %')
  return(results)
}

ndw_lasso_over_p <- function(n, list_p, s_t, delta, ...){
  n_iterations = length(list_p)
  results = matrix(rep(0, n_iterations*6), nrow = n_iterations)
  
  for(i in 1:n_iterations){
    corr_matrix = get_correlation_matrix(list_p[i], s_t, delta)
    true_edge_set = corr_matrix != 0
    sample = generate_multivariate_norm_sample(n, corr_matrix)
    edge_set = ndw_lasso(sample, ...)$edge_set
    results[i,] = summary_stat(true_edge_set, edge_set)
  }
  colnames(results) = c('TPR', 'FPR', 'TNR', 'FNR', 'TC %', 'Error %')
  return(results)
}

ndw_lasso_over_s_t <- function(n, p, list_s_t, delta, ...){
  n_iterations = length(list_s_t)
  results = matrix(rep(0, n_iterations*6), nrow = n_iterations)
  
  for(i in 1:n_iterations){
    corr_matrix = get_correlation_matrix(p, list_s_t[i], delta)
    true_edge_set = corr_matrix != 0
    sample = generate_multivariate_norm_sample(n, corr_matrix)
    edge_set = ndw_lasso(sample, ...)$edge_set
    results[i,] = summary_stat(true_edge_set, edge_set)
  }
  colnames(results) = c('TPR', 'FPR', 'TNR', 'FNR', 'TC %', 'Error %')
  return(results)
}

# GLASSO

# Glasso is computed directly using the library glasso usually, so here we use that 
# and the wrapper Cvglasso instead of creating an explicit glasso function.

# Repeats the glasso 50 times, selecting optimal regularization using k-fold
# cross validation. 
# K controls number of folds
# If cross_validate = FALSE, constructs a set of lambdas and gives edge set
# using the best possible lambda. 
# crit.cv controls which validation method to use (loglikelihood, AIC, BIC)
# nlambda controls number of lambdas to cross-validate with
# scale controls whether or not data is standardized prior to applying glasso.
# ... <- additional arguments to pass to CV glasso.
# Check CVglasso documentation for more info

glasso_repeater <- function(n, p , s_t, delta, times_repeated= 50, K = 5,
                            cross_validate=TRUE, crit.cv = 'loglik', nlambda = 30, 
                            lambda_min = 0.01, lambda_max = 1,lam.min.ratio = 0.0001,
                            scale = FALSE,  ...){
  results = matrix(rep(0, times_repeated * 6), nrow = times_repeated)
  sparsity_pattern = get_correlation_matrix(p, s_t, delta)
  lambdas = seq(from = lambda_min, to = lambda_max, length.out = nlambda)
  true_edge_set = sparsity_pattern != 0
  progress_bar = txtProgressBar(min = 0, max = times_repeated, style =3)
  if(cross_validate){
    for (repitition in 1:times_repeated){
      generated_sample = generate_multivariate_norm_sample(n, sparsity_pattern)
      if(scale){
        generated_sample = scale(generated_sample)
      }
      cross_validated_result = CVglasso(X=generated_sample, nlam=nlambda, 
                                        lam.min.ratio = lam.min.ratio, K=K,
                                        crit.cv = crit.cv, trace="none", ...)$Omega
      edge_set = cross_validated_result != 0
      summary = summary_stat(true_edge_set, edge_set, diagonal=TRUE)
      results[repitition, ] = summary
      setTxtProgressBar(progress_bar, repitition)
    }
  }
  else{
    for(repitition in 1:times_repeated){
      generated_sample = generate_multivariate_norm_sample(n, sparsity_pattern)
      if(scale){
        generated_sample= scale(generated_sample)
      }
      covar_sample = cov(generated_sample)
      glasso_results = glassopath(covar_sample, rholist = lambdas, 
                                  trace = 0, ...)$wi!=0
      summary = apply(glasso_results, MARGIN = 3, FUN = summary_stat, 
                      true_edge_set=true_edge_set, diagonal = TRUE)
      results[repitition,] = summary[ ,which.min(summary[6,])]
      setTxtProgressBar(progress_bar, repitition)
    }
  }

  close(progress_bar)
  colnames(results) = c('TPR', 'FPR', 'TNR', 'FNR', 'TC %', 'Error %')
  return(results)
}

glasso_over_n <-function(list_n, p, s_t, delta, ...) {
  n_iterations = length(list_n)
  results = matrix(rep(0, n_iterations*6), nrow = n_iterations)
  corr_matrix = get_correlation_matrix(p, s_t, delta)
  true_edge_set = corr_matrix != 0
  for(i in 1:n_iterations){
    sample = generate_multivariate_norm_sample(list_n[i], corr_matrix)
    glasso_results = CVglasso(X=sample, trace='none', ...)$Omega
    edge_set = glasso_results != 0
    results[i,] = summary_stat(true_edge_set, edge_set, diagonal = TRUE)
  }
  colnames(results) = c('TPR', 'FPR', 'TNR', 'FNR', 'TC %', 'Error %')
  return(results)
}

glasso_over_p <- function(n, list_p, s_t, delta, ...){ 
  n_iterations = length(list_p)
  results = matrix(rep(0, n_iterations*6), nrow = n_iterations)
  progress_bar = txtProgressBar(min = 0, max = n_iterations, style =3)
  
  for(i in 1:n_iterations){
    corr_matrix = get_correlation_matrix(list_p[i], s_t, delta)
    true_edge_set = corr_matrix != 0
    sample = generate_multivariate_norm_sample(n, corr_matrix)
    glasso_results = CVglasso(X=sample, trace='none',...)$Omega
    edge_set = glasso_results != 0
    results[i,] = summary_stat(true_edge_set, edge_set, diagonal = TRUE)
    setTxtProgressBar(progress_bar, i)
  }
  close(progress_bar)
  colnames(results) = c('TPR', 'FPR', 'TNR', 'FNR', 'TC %', 'Error %')
  return(results)
}

glasso_over_s_t <- function(n, p, list_s_t, delta, ...){ # s_t controls sparsity structure
  n_iterations = length(list_s_t)
  results = matrix(rep(0, n_iterations*6), nrow = n_iterations)
  
  for(i in 1:n_iterations){
    corr_matrix = get_correlation_matrix(p, list_s_t[i], delta)
    true_edge_set = corr_matrix != 0
    sample = generate_multivariate_norm_sample(n, corr_matrix)
    glasso_results = CVglasso(X=sample, trace='none', ...)$Omega
    edge_set = glasso_results != 0
    results[i,] = summary_stat(true_edge_set, edge_set, diagonal = TRUE)
  }
  colnames(results) = c('TPR', 'FPR', 'TNR', 'FNR', 'TC %', 'Error %')
  return(results)
}
# END INTIALIZATION

# 1.GETTING ROC CURVES
# Change values and re run with same seed to replicate report results
# 1.1 Make Sample and get a list of lambdas

n = 500 #100 #50
p = 100
s_t = 0.9 #0.6
delta = 2.01

set.seed(101)
data = lasso_sample_maker(n, p, s_t, delta)
sample = data$sample
true_edge_set = data$sparsity_pattern != 0

lambda_ndw = seq(from = 5, to = 0.001, length.out = 5000) #ndw requires decreasing sequence
lambda_glasso = rev(lambda_ndw) #glassopath requires increasing sequence

# 1.2 Nodewise Lasso AUROC, using two different kinds of edge retrieval 
edge_sets_and = ndw_lasso(sample, lambdalist = lambda_ndw, 
                          method = 'lambdalist', and = TRUE)$edge_set

edge_sets_or = ndw_lasso(sample, lambdalist = lambda_ndw, 
                         method = 'lambdalist', and = FALSE)$edge_set 


summary_stats_and = t(sapply(edge_sets_and, summary_stat, 
                           true_edge_set = true_edge_set))

summary_stats_or = t(sapply(edge_sets_or, summary_stat, 
                          true_edge_set = true_edge_set))

AUC_and  = signif(auc(summary_stats_and[,2], summary_stats_and[,1], 
                      type = 'spline'), 4)
AUC_or  = signif(auc(summary_stats_or[,2], summary_stats_or[,1], 
                     type = 'spline'), 4)

par(mfrow = c(1,2))

plot(summary_stats_and[,2], summary_stats_and[,1], ylab= 'TPR', xlab = 'FPR',
     main = 'ROC curve (and)', type = 'l', xlim=c(0, 1), ylim = c(0,1))
abline(0,1)
text(0.5, 0.0, paste('Area under curve:', AUC_and))

plot(summary_stats_or[,2], summary_stats_or[,1], ylab= 'TPR', xlab = 'FPR', 
     main = 'ROC curve (or)', type = 'l', xlim= c(0,1), ylim = c(0,1))
abline(0,1)
text(0.5, 0.0, paste('Area under curve:', AUC_or))

# 1.3 Glasso AUROC
set.seed(103)

par(mfrow = c(1,1))

edge_sets_glasso = glassopath(cov(sample), rholist = lambda_glasso, trace = 0)$wi != 0

summary_glasso = t(apply(edge_sets_glasso, MARGIN = 3, FUN = summary_stat, 
                true_edge_set=true_edge_set, diagonal = TRUE))

AUC_glasso = signif(auc(summary_glasso[,2], summary_glasso[,1], 
                        type = 'spline'), 4)

plot(summary_glasso[,2], summary_glasso[,1], ylab= 'TPR', xlab = 'FPR',
     main = 'ROC curve (Glasso)', type = 'l', xlim=c(0, 1), ylim = c(0,1))
abline(0,1)
text(0.5, 0.0, paste('Area under curve:', AUC_glasso))

# 1.4 Change in number of predicted edges with lambda for glasso

edge_pred = apply(edge_sets_glasso, 3, function(e) (sum(1*(e==TRUE))-p)/2)
plot(lambda_glasso, edge_pred, ylab = 'Number of Edges Predicted', type = 'l', main='Glasso')

# 1.5 Change in prediction error with lambdas:

plot(lambda_ndw,summary_stats_or[,6], type ='l', ylab= 'Error', main= 'Node-wise Lasso')

plot(lambda_glasso, summary_glasso[,6], type ='l', ylab= 'Error', main='Glasso')

# 2. Repeating 50 times
#  2.1 Comparing MBpen and CV for Node-wise Lasso
# Change values and re run with same seed to replicate report results

n = 500 #50 #100
p = 100
s_t = 0.9
delta = 1

set.seed(201)

# 2.1.1 Nodewise Lasso (AND) with Meinshausen-Buhlmann Optimal Lambda
par(mfrow = c(1,2))

results_ndw_mb = ndw_lasso_repeater(n,p, s_t, delta, and=TRUE, method ='mbpen', alpha =0.05)
boxplot(results_ndw_mb, ylab = 'Rate/Percent', main = 'Node-wise Lasso with MB Penalty')
colMeans(results_ndw_mb) # Mean Values
sqrt(diag(var(results_ndw_mb))/50) #Standard Errors


# 2.1.2 Nodewise Lasso (AND) with 5-fold cross validation

results_ndw_cv = ndw_lasso_repeater(n,p,s_t,delta, and = TRUE, method = 'cv', k = 5)
boxplot(results_ndw_cv, ylab = 'Rate/Percent', main = 'Node-wise Lasso with 5-fold CV')
colMeans(results_ndw_cv) # Mean Values
sqrt(diag(var(results_ndw_cv))/50) #Standard Errors



# 2.2 Results with Optimal Parameter Selection, MBpen versus Graphical Lasso (Repeated 50 times)
# Change values and re run with same seed to replicate report results
n = 500 #100
p = 100
s_t = 0.9 #0.6
delta = 1
par(mfrow = c(1,2))

set.seed(500)

# 2.2.1 Nodewise Lasso (AND) with Meinshausen-Buhlmann Optimal Lambda

results_ndw_mb2 = ndw_lasso_repeater(n,p, s_t, delta, and=TRUE, method ='mbpen', alpha =0.05)
boxplot(results_ndw_mb2, ylab = 'Rate/Percent', main = 'Node-wise Lasso with MB Penalty')
colMeans(results_ndw_mb2) # Mean Values
sqrt(diag(var(results_ndw_mb2))/50) #Standard Errors


# 2.2.2 Glasso with 5-fold cross validation 
# Slow 
set.seed(500)
results_glas_5fold = glasso_repeater(n,p, s_t, delta, nlambda= 30)
boxplot(results_glas_5fold, ylab = 'Rate/Percent', main = 'Glasso with 5-fold CV')
sqrt(diag(var(results_glas_5fold))/50) #Standard Errors
colMeans(results_glas_5fold) # Mean Values


# 3. Varying n, p, s_t
# 3.1 Results with variation in n (Sparse Matrix):

p = 50
s_t = 0.9
delta = 1
list_n = seq(from = 50, to = 500, by = 20)
set.seed(950)
results_ndw_n = ndw_lasso_over_n(list_n,p,s_t,delta)
results_glasso_n = glasso_over_n(list_n,p,s_t,delta, lam.min.ratio = 0.0001)

par(mfrow = (c(1,2)))
plot(list_n, results_ndw_n[,6], ylab = 'Error %/ True Positive Rate',
     main = 'Node-wise Lasso over n', col='red', ylim = c(0,1), xlab = 'n', pch=20)
points(list_n, results_ndw_n[,1], col='green', pch=20)
legend(310,0.8, c('Error %', 'TPR'), col = c('red', 'green'), pch = c(20,20))

plot(list_n, results_glasso_n[,6], ylab = 'Error %/ True Positive Rate',
     main = 'Glasso over n', col = 'red', ylim = c(0,1), xlab= 'n', pch=20)
points(list_n, results_glasso_n[,1], col='green', pch=20)
legend(310,0.8, c('Error %', 'TPR'), col = c('red', 'green'), pch = c(20,20))



# 3.2 Results with variation in n (Non-sparse Matrix):

p = 50
s_t = 0.6 # Very dense matrix, P[conditionally independent] = 0.3.
delta = 1
list_n = seq(from = 50, to = 500, by = 20)
set.seed(950)

results_ndw_n = ndw_lasso_over_n(list_n,p,s_t,delta)
results_glasso_n = glasso_over_n(list_n,p,s_t,delta, lam.min.ratio = 0.0001)

par(mfrow = (c(1,2)))
plot(list_n, results_ndw_n[,6], ylab = 'Error %/ True Positive Rate',
     main = 'Node-wise Lasso over n', col='red', ylim = c(0,1), xlab = 'n', pch=20)
points(list_n, results_ndw_n[,1], col='green', pch=20)
legend(310,0.8, c('Error %', 'TPR'), col = c('red', 'green'), pch = c(20,20))

plot(list_n, results_glasso_n[,6], ylab = 'Error %/ True Positive Rate',
     main = 'Glasso over n', col = 'red', ylim = c(0,1), xlab= 'n', pch=20)
points(list_n, results_glasso_n[,1], col='green', pch=20)
legend(310,0.8, c('Error %', 'TPR'), col = c('red', 'green'), pch = c(20,20))


# 3.3 Results with Variation in p (Sparse Matrix):

n = 500
s_t = 0.9
delta = 2.01
list_p = seq(from= 10, to= 230, by = 20)

set.seed(100)
results_ndw_p = ndw_lasso_over_p(n,list_p,s_t,delta)
results_glasso_p = glasso_over_p(n,list_p,s_t,delta, lam.min.ratio = 0.0001) 

par(mfrow = (c(1,2)))
plot(list_p, results_ndw_p[,6], ylab = 'Error %/ True Positive Rate',
     main = 'Node-wise Lasso over p', col='red', ylim = c(0,1), xlab = 'p', pch=20)
points(list_p, results_ndw_p[,1], col='green', pch=20)
legend(160,0.8, c('Error %', 'TPR'), col = c('red', 'green'), pch = c(20,20))

plot(list_p, results_glasso_p[,6], ylab = 'Error %/ True Positive Rate',
     main = 'Glasso over p', col = 'red', ylim = c(0,1), xlab= 'p', pch=20)
points(list_p, results_glasso_p[,1], col='green', pch=20)
legend(160,0.8, c('Error %', 'TPR'), col = c('red', 'green'), pch = c(20,20))


# 3.4 Results with variation in sparsity (given n, p):

n = 1000
p = 50
delta = 1
list_s_t = seq(from = 0.1, to = 0.9, by = 0.1)

set.seed(150)
results_ndw_s_t_mbpen = ndw_lasso_over_s_t(n,p,list_s_t,delta)
results_ndw_s_t_cv = ndw_lasso_over_s_t(n,p,list_s_t,delta, method = 'cv')
results_glasso_s_t = glasso_over_s_t(n,p,list_s_t,delta, lam.min.ratio = 0.0001)
# Will say tuning parameter on boundary, but that is natural because in dense settings
# an optimal prediction is when lambda is nearly 0.

par(mfrow = (c(1,2)))

plot(list_s_t, results_ndw_s_t_mbpen[,6], ylab = 'Error %/ True Positive Rate',
     main = 'Node-wise Lasso over s_t (MBpen)', col='red', ylim = c(0,1),
     xlab = 'Sparsity', pch=20)
points(list_s_t, results_ndw_s_t_mbpen[,1], col='green', pch=20)
legend(0.6,0.8, c('Error %', 'TPR'), col = c('red', 'green'), pch = c(20,20))

plot(list_s_t, results_ndw_s_t_cv[,6], ylab = 'Error %/ True Positive Rate',
     main = 'Node-wise Lasso over s_t (CV)', col='red', ylim = c(0,1), 
     xlab = 'Sparsity', pch=20)
points(list_s_t, results_ndw_s_t_cv[,1], col='green', pch=20)
legend(0.6,0.8, c('Error %', 'TPR'), col = c('red', 'green'), pch = c(20,20))

par(mfrow = c(1,1))

plot(list_s_t, results_glasso_s_t[,6], ylab = 'Error %/ True Positive Rate',
     main = 'Glasso over s_t', col = 'red', ylim = c(0,1), xlab= 'Sparsity', pch=20)
points(list_s_t, results_glasso_s_t[,1], col='green', pch=20)
legend(0.6,0.8, c('Error %', 'TPR'), col = c('red', 'green'), pch = c(20,20))


# 4. Variation in data generating process (Adding uniform noise to data): 
# Run block to replicate results
n = 500
p = 100
s_t = 0.9
delta = 2.01

set.seed(100)
data = lasso_sample_maker(n, p, s_t,delta)
Z = data$sample
true_edge_set = data$sparsity_pattern != 0
Y = matrix(runif(n*p,-1,1), nrow = n, ncol = p)
X = Z + Y

edge_set_ndw_z = ndw_lasso(Z)$edge_set
edge_set_glasso_z = CVglasso(Z, lam.min.ratio = 0.0001, trace = 'none')$Omega != 0

edge_set_ndw_x = ndw_lasso(X)$edge_set
edge_set_glasso_x = CVglasso(X, lam.min.ratio = 0.0001, trace = 'none')$Omega != 0


ndw_results_z = summary_stat(true_edge_set, edge_set_ndw_z)
glasso_results_z = summary_stat(true_edge_set, edge_set_glasso_z, diagonal = TRUE)

ndw_results_x = summary_stat(true_edge_set, edge_set_ndw_x)
glasso_results_x = summary_stat(true_edge_set, edge_set_glasso_x, diagonal = TRUE)

ndw_results_z
ndw_results_x

glasso_results_z
glasso_results_x

# End Of Code
# Thanks for reading

