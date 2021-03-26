# Estimator 3: CAUSAL FOREST
# Reference:  Wager and Athey [2018]; Athey et al. [2019]
#

CATE_CausalForest <- function(Y,W,X,X_test){
  
  # STEP 1: Fit the forest
  ## Fit the forest (Double-sample Trees procedure)
  cf <- grf::causal_forest(X = as.matrix(X),
                           Y = as.matrix(Y),
                           W = as.matrix(W),
                           num.trees = 2000, # Make this larger for better acc.
                           honesty = TRUE,
                           min.node.size = 50)
  
  # STEP 2.a: Predict point estimates and standard errors (training set, out-of-bag)
  #oob_pred <- predict(cf, estimate.variance=TRUE)
  #oob_tauhat_cf <- oob_pred$predictions
  #oob_tauhat_cf_se <- sqrt(oob_pred$variance.estimates)
  
  # STEP 2.b: Predict point estimates and standard errors (test set)
  test_pred <- predict(cf, newdata = X_test, estimate.variance=TRUE)
  tauhat_cf <- test_pred$predictions
  #sehat_cf <- sqrt(test_pred$variance.estimates)
  
  return(tauhat_cf)
}