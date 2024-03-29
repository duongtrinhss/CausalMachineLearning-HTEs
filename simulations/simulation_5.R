# SIMULATION 5: Low dimensional data; Confounding;
# Linear, sparse mean effect; Nonlinear treatment effect
#

simulation_5 <- function(d){
  # Generate training sample
  n <- 5000
  sigma <- 1
  
  ## Covariates X
  X <-  matrix(runif(n * d, 0, 1), n, d)  # features
  
  ## Propensity function: 
  propensity <- function(x){
    save_the_seed_1 <- .Random.seed
    set.seed(208)
    cf <- dbeta(x[1], 2, 4)
    .Random.seed <<- save_the_seed_1
    
    return(((1 + cf)/4))
    
    }
  W <-  rbinom(n, 1, apply(X,1,propensity))
  
  ## Mean effect function: 
  mean_effect <- function(x){
    
    (2*x[1] -1)
    
  }
  
  ## Treatment effect function: tau(x)
  treatment_effect <-  function(x) {
    
    (1 + 1 / (1 + exp(-20 * (x[1] - 1/3)))) * (1 + 1 / (1 + exp(-20 * (x[2] - 1/3))))
    
  } # Equation (28) for heterogeneity effect
  
  ## Outcome
  Y <-  (W - 0.5) * apply(X, 1, treatment_effect) + sigma * rnorm(n)  
  
  # Generate testing sample
  n_test <- 1000
  save_the_seed_2 <- .Random.seed
  set.seed(261) # To maintain the testing sample through simulations
  X_test <-  matrix(runif(n_test * d, 0, 1), n_test, d)
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
  print("DRLearner")
  tauhat_dr <- CATE_DRLearner(Y,W,X,X_test)
  
  tau_all <- data.frame(true_effect, tauhat_pl,tauhat_ct,tauhat_cf,tauhat_xl,tauhat_dr)
  
  return(tau_all)
}

# Simulation for Causal Forest local centering

simulation_5b <- function(d){
  # Generate training sample
  n <- 5000
  sigma <- 1
  
  ## Covariates X
  X <-  matrix(runif(n * d, 0, 1), n, d)  # features
  
  ## Propensity function: 
  propensity <- function(x){
    save_the_seed_1 <- .Random.seed
    set.seed(208)
    cf <- dbeta(x[1], 2, 4)
    .Random.seed <<- save_the_seed_1
    
    return(((1 + cf)/4))
    
  }
  W <-  rbinom(n, 1, apply(X,1,propensity))
  
  ## Mean effect function: 
  mean_effect <- function(x){
    
    (2*x[1] -1)
    
  }
  
  ## Treatment effect function: tau(x)
  treatment_effect <-  function(x) {
    
    (1 + 1 / (1 + exp(-20 * (x[1] - 1/3)))) * (1 + 1 / (1 + exp(-20 * (x[2] - 1/3))))
    
  } # Equation (28) for heterogeneity effect
  
  ## Outcome
  Y <-  (W - 0.5) * apply(X, 1, treatment_effect) + sigma * rnorm(n)  
  
  # Generate testing sample
  n_test <- 1000
  save_the_seed_2 <- .Random.seed
  set.seed(261) # To maintain the testing sample through simulations
  X_test <-  matrix(runif(n_test * d, 0, 1), n_test, d)
  .Random.seed <<- save_the_seed_2
  X_test <- as.data.frame(X_test)
  ## Calculate true CATE on testing data for evaluation purpose (calculating expected MSE)
  true_effect <-  apply(X_test, 1, treatment_effect)
  
  # Predict CATEs
  print("CausalForest local centering")
  tauhat_cf_lc <- CATE_CausalForest_lc(Y,W,X,X_test)

  
  tau_all <- data.frame(true_effect, tauhat_cf_lc)
  
  return(tau_all)
}
