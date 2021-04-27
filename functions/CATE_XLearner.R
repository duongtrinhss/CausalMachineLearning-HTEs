# Estimator 4: X - LEARNER
# Reference:  Kuenzel [2019]
#

CATE_XLearner <- function(Y,W,X,X_test){
  
  num.trees <- 2000
  
  tauhat_xl = matrix(0, nrow(X_test), 1)
  
  # Create folds for Cross-fitting
  index <-  caret::createFolds(Y, k = 2)
  
  for(i in 1:length(index)) {
    
    # STEP 1: Fit the X-learner
  
    ## Estimate RF nuisance parameters
    ### Estimate separate response functions
    tf0 <- regression_forest(X[-index[[i]],,drop=F][W[-index[[i]]]==0,], Y[-index[[i]]][W[-index[[i]]]==0], num.trees = num.trees) # RF base learner for the first stage
    tf1 <- regression_forest(X[-index[[i]],,drop=F][W[-index[[i]]]==1,], Y[-index[[i]]][W[-index[[i]]]==1], num.trees = num.trees)
    ### Estimate the propensity score
    propf <- regression_forest(X[-index[[i]],,drop=F], W[-index[[i]]], num.trees = num.trees)
  
  
    ## Compute the "imputed treatment effects" using the other group
    W0 <- predict(tf1, X[index[[i]],,drop=F][W[index[[i]]]==0,])$predictions - Y[index[[i]]][W[index[[i]]]==0]
    W1 <- Y[index[[i]]][W[index[[i]]]==1] - predict(tf0, X[index[[i]],,drop=F][W[index[[i]]]==1,])$predictions
  
   ## Compute the cross estimators
    xf0 <- regression_forest(X[index[[i]],,drop=F][W[index[[i]]]==0,], W0, num.trees = num.trees)
    xf1 <- regression_forest(X[index[[i]],,drop=F][W[index[[i]]]==1,], W1, num.trees = num.trees)
  
  
   # STEP 2: Predict point estimates on the test set
    ehat.test <- predict(propf, X_test)$predictions
    xf.preds.0.test <- predict(xf0, X_test)$predictions
    xf.preds.1.test <- predict(xf1, X_test)$predictions
    tauhat_idx <- (1-ehat.test)*xf.preds.1.test + ehat.test*xf.preds.0.test
    
    tauhat_xl <- tauhat_xl + 1/length(index)*tauhat_idx
    
  }
  
  return(tauhat_xl)
  
  
  
}

# CATE_XLearner <- function(Y,W,X,X_test){
#   
#   num.trees <- 2000  #  We'll make this a small number for speed here.
#   
#   # STEP 1: Fit the X-learner
#   ## Estimate separate response functions
#   tf0 <- regression_forest(X[W==0,], Y[W==0], num.trees = num.trees) # RF base learner for the first stage
#   tf1 <- regression_forest(X[W==1,], Y[W==1], num.trees = num.trees)
#   ## Compute the "imputed treatment effects" using the other group
#   D1 <- Y[W==1] - predict(tf0, X[W==1,])$predictions
#   D0 <- predict(tf1, X[W==0,])$predictions - Y[W==0]
#   ## Compute the cross estimators
#   xf0 <- regression_forest(X[W==0,], D0, num.trees = num.trees)
#   xf1 <- regression_forest(X[W==1,], D1, num.trees = num.trees)
#   ## Predict treatment effects, making sure to always use OOB predictions where appropriate
#   xf.preds.0 <- rep(0,length(Y))
#   xf.preds.0[W==0] <- predict(xf0)$predictions # ? why not predict(xf0, X[W==0,])$predictions - out of bag?
#   xf.preds.0[W==1] <- predict(xf0, X[W==1,])$predictions 
#   
#   xf.preds.1 <- rep(0,length(Y))
#   xf.preds.1[W==0] <- predict(xf1, X[W==0,])$predictions 
#   xf.preds.1[W==1] <- predict(xf1)$predictions 
#   ## Estimate the propensity score
#   propf <- regression_forest(X, W, num.trees = num.trees)
#   ehat <- predict(propf)$predictions
#   ## Finally, compute the X-learner prediction
#   oob_tauhat_xl <- (1-ehat)*xf.preds.1 + ehat*xf.preds.0
#   
#   # STEP 2: Predict point estimates on the test set
#   ehat.test <- predict(propf, X_test)$predictions
#   xf.preds.0.test <- predict(xf0, X_test)$predictions
#   xf.preds.1.test <- predict(xf1, X_test)$predictions
#   tauhat_xl <- (1-ehat.test)*xf.preds.1.test + ehat.test*xf.preds.0.test
#   
#   return(tauhat_xl)
# }

# STEP 3: Compute standard error - Bootstrap for Normal approximated CI:
#   computeSE <- function(Y,W,X,X_test,B){
#     s0 <- seq(length(W))[W==0]
#     s1 <- seq(length(W))[W==1]
#     n0 <- length(s0) 
#     n1 <- length(s1)
#     Tauhat_xl_test <- matrix(rep(0,dim(X_test)[1]*B),nrow = dim(X_test)[1])
#     for(b in 1:B){
#       sb <- c(sample(s0, n0, replace = TRUE),sample(s1, n1, replace = TRUE))
#       yb = Y[sb]
#       xb = X[sb,]
#       wb = W[sb]
#       Tauhat_xl_test[,b] <- X_learner(yb,wb,xb,X_test)
#     } 
#     return(apply(Tauhat_xl_test,1,sd))
#   }
# sehat_xl <- computeSE(Y,W,X,X_test,10)