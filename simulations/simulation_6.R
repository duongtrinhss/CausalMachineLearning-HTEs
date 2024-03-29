# SIMULATION 6: High dimensional data; No Confounding; Balanced propensity propensity; 
# No mean effect; Nonlinear treatment effect
#

simulation_6 <- function(d){
  # Generate training sample
  n <- 200
  sigma <- 1
  
  ## Covariates X
  X <-  matrix(runif(n * d, 0, 1), n, d)  # features
  
  ## Propensity function: e(x) = 0.5 - constant
  W <-  rbinom(n, 1, 0.5)
  
  ## Mean effect function: m(x) = 0 - constant
  
  ## Treatment effect function: tau(x)
  treatment_effect <-  function(x) {
    (1 + 1 / (1 + exp(-20 * (x[1] - 1/3)))) * (1 + 1 / (1 + exp(-20 * (x[2] - 1/3))))
  } # Equation (28) for heterogeneity effect
  
  Y <-  (W - 0.5) * apply(X, 1, treatment_effect) + sigma * rnorm(n)  
  
  # Generate testing sample
  n_test <- 200
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
