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


CATE_CausalForest_lc <- function(Y,W,X,X_test){
  
  num.trees <- 2000
  
  tauhat_cf = matrix(0, nrow(X_test), 1)
  np <-  matrix(NA,length(W),2)
  
  # Create folds 
  index <-  caret::createFolds(Y, k = 2)
  
  for(i in 1:length(index)) {
    # P-score
    fit_p <- do.call(regression_forest, c(list(X=X[-index[[i]],,drop=F],
                                               Y=W[-index[[i]]]),
                                          #tune.parameters = "none",
                                          num.trees = num.trees))
    np[index[[i]],1] <-  predict(fit_p,X[index[[i]],,drop=F])$prediction
    
    # Outcome
    fit_y <- do.call(regression_forest, c(list(X=X[-index[[i]],,drop=F],
                                               Y=Y[-index[[i]]]),
                                          #tune.parameters = "none",
                                          num.trees = num.trees))
    np[index[[i]],2] <-  predict(fit_y,X[index[[i]],,drop=F])$prediction
  }
    
    
    # STEP 1: Fit the forest
    ## Fit the forest (Double-sample Trees procedure)
    cf <- grf::causal_forest(X = as.matrix(X),
                             Y = as.matrix(Y),
                             W = as.matrix(W),
                             W.hat = np[,1],
                             Y.hat = np[,2],
                             num.trees = 2000, # Make this larger for better acc.
                             honesty = TRUE,
                             min.node.size = 50)
    
    # STEP 2.b: Predict point estimates and standard errors (test set)
    test_pred <- predict(cf, newdata = X_test, estimate.variance=TRUE)
    tauhat_cf <- test_pred$predictions
    
  return(tauhat_cf)
}
