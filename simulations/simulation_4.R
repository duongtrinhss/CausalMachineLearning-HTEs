# SIMULATION 4: Low dimensional data; No confounding; Unbalanced propensity propensity; 
# Linear, dense/sparse mean effect; Linear treatment effect
#

simulation_4 <- function(d){
  # Generate training sample
  n <- 5000
  sigma <- 1
  
  ## Covariates X
  ## Draw a random nxd matrix from MN(M,U,V)
  Mean_mat <- matrix(0, nrow = n, ncol = d)
  U <- I(n)
  V <- I(d)
  #X <- rmatnorm(s=1,Mean.mat,U,V)  # features
  X <- rmatrixnorm(n=1,mean=Mean_mat,U=U,V=V)
  ## Propensity function: e(x) = 0.5 - constant
  W <-  rbinom(n, 1, 0.01)
  
  ## Mean effect function: 
  mean_effect <- function(x){
    save_the_seed_1 <- .Random.seed
    set.seed(208) # To maintain through simulations   
    b <- runif(d, -5, 5)
    .Random.seed <<- save_the_seed_1
    
    return(x%*%b)
    
  }
  
  ## Treatment effect function: tau(x)
  treatment_effect <-  function(x){
    (6*(x[1]>0)*x[1] + 8*(x[2]>0))*x[2]
  }
  Y <- apply(X, 1, mean_effect) + (W - 0.5) * apply(X, 1, treatment_effect) + sigma * rnorm(n)  
  
  # Generate testing sample
  n_test <- 1000
  Mean_mat_test <- matrix(0, nrow = n_test, ncol = d)
  U_test <- I(n_test)
  V_test <- I(d)
  #X_test <- rmatnorm(s=1,Mean.mat.test,U.test,V.test)
  save_the_seed_2 <- .Random.seed
  set.seed(261) # To maintain the testing sample through simulations
  X_test <- rmatrixnorm(n=1,mean=Mean_mat_test,U=U_test,V=V_test)
  .Random.seed <<- save_the_seed_2
  X_test <- as.data.frame(X_test)
  
  ## Calculate true CATE on testing data for evaluation purpose (calculating expected MSE)
  true_effect <-  apply(X_test, 1, treatment_effect)
  
  # Predict CATEs
  print("PostLasso")
  tauhat_pl <- CATE_PostLasso(Y,W,X,X_test)
  print("CausalTree")
  tauhat_ct <- CATE_CausalTree(Y,W,X,X_test)
  print("CausalForest")
  tauhat_cf <- CATE_CausalForest(Y,W,X,X_test)
  print("XLearner")
  tauhat_xl <- CATE_XLearner(Y,W,X,X_test)
  print("MomGrf")
  tauhat_mom <- CATE_MomGrf(Y,W,X,X_test)
  
  tau_all <- data.frame(true_effect, tauhat_pl,tauhat_ct,tauhat_cf,tauhat_xl,tauhat_mom)
  
  return(tau_all)
}
