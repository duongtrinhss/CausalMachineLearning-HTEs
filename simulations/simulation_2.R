# SIMULATION 2: Low dimensional data; No confounding; Balanced propensity; 
# Nonlinear mean effect; Linear treatment effect
#

simulation_2 <- function(d){
  # Generate training sample
  n <- 5000
  sigma <- 1
  
  ## Covariates X
  X <-  matrix(runif(n * d, 0, 1), n, d)  # features
  
  ## Propensity function: e(x) = 0.5 - constant
  W <-  rbinom(n, 1, 0.5)
  
  ## Mean effect function: 
  mean_effect <- function(x){
    sin(pi*x[1]*x[2]) + 2*(x[3] - 0.5)^2 + x[4] + 0.5*x[5]
  }
  
  ## Treatment effect function: tau(x)
  treatment_effect <-  function(x) {
    1/2*(x[1]+x[2])
  } 
    Y <- apply(X, 1, mean_effect) + (W - 0.5) * apply(X, 1, treatment_effect) + sigma * rnorm(n)  
  
  # Generate testing sample
  n_test <- 1000
  save_the_seed <- .Random.seed
  set.seed(261) # To maintain the testing sample through simulations
  X_test <-  matrix(runif(n_test * d, 0, 1), n_test, d)
  .Random.seed <<- save_the_seed
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
  print("DRLearner")
  tauhat_dr <- CATE_DRLearner(Y,W,X,X_test)
  
  tau_all <- data.frame(true_effect, tauhat_pl,tauhat_ct,tauhat_cf,tauhat_xl,tauhat_dr)
  
  return(tau_all)
}