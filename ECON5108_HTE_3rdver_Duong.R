#===============================================================================
# Written on March, 2021
# Duong Trinh
# ECON5108: Issues in Economic Research
# MRes in Economics - University of Glasgow
# Estimation of Heterogeneous Treatment Effects  
#===============================================================================

# setwd("C/Users/T450s/Documents/projects/ECON5108-DT")

rm(list = ls())
# Loading packages
source("packages.R")

#===============================| NOTATION NOTES |==============================

# Y, observed outcome on training sample
# W, binary indicator for the treatment on training sample
# X, observed variables or feature matrix on training sample
# X.test, observed variables on testing sample

#====================| Causal Machine Learning Estimators  |====================

# Estimator 1: POST LASSO
#
source("functions/CATE_PostLasso.R")

# Estimator 2: (HONEST) CAUSAL TREE 
# Reference:  Athey and Imbens [2016]
#
source("functions/CATE_CausalTree.R")

# Estimator 3: CAUSAL FOREST
# Reference:  Wager and Athey [2018]; Athey et al. [2019]
#
source("functions/CATE_CausalForest.R")

# Estimator 4: X - LEARNER
# Reference:  Kuenzel [2019]
#
source("functions/CATE_XLearner.R")

# Estimator 5: MOM - GRF
# Reference:  Knaus et al. [2021]
#
source("functions/CATE_MomGrf.R")

#==============| COMPARISON METRICS FOR 4 CAUSAL ML ESTIMATORS |================

#=============================| SIMULATION DESIGN |=============================

# SIMULATION SETUP 1: Low dimensional data; No confounding; 
# Balanced propensity; Linear mean effect; No treatment effect
#
source("simulations/simulation_1.R")

# SIMULATION 2: Low dimensional data; No confounding; Balanced propensity; 
# Linear, dense/sparse mean effect; Linear treatment effect
#
source("simulations/simulation_2.R")

# SIMULATION 3: Low dimensional data; No confounding; Balanced propensity; 
# No mean effect; Nonlinear treatment effect
#
# Reference: Wager & Athey(2017)
#
source("simulations/simulation_3.R")

# SIMULATION 4: Low dimensional data; No confounding; Unbalanced propensity propensity; 
# Linear, dense/sparse mean effect; Linear treatment effect
#
source("simulations/simulation_4.R")

# SIMULATION 5: Low dimensional data; Confounding; Balanced propensity propensity; 
# Linear, sparse mean effect; Nonlinear treatment effect
#
source("simulations/simulation_5.R")

# SIMULATION 6: High dimensional data; No Confounding; Balanced propensity propensity; 
# No mean effect; Nonlinear treatment effect
#
source("simulations/simulation_6.R")



#====================| COLLECT RESULTS OF MC SIMULATIONS|=======================
# Parameters for Simulation 1
dvals <- c(5)
iter <- 30
n_test <- 1000

# Parameters for Simulation 2
dvals <- c(10,20)
iter <- 25
n_test <- 1000

# Parameters for Simulation 3
dvals <- c(2,4,6,8)
iter <- 25
n_test <- 1000

# Parameters for Simulation 4
dvals <- c(2,8)
iter <- 25
n_test <- 1000

# Parameters for Simulation 5
dvals <- c(2,8)
iter <- 25
n_test <- 1000

# Parameters for Simulation 6
dvals <- c(200)
iter <- 25
n_test <- 200

# Run simulations in turn by replacing line 105
result_raw <- lapply(dvals, function(d){
  print(paste("NOW RUNNING: d = ", d))
  # Get CATEs matrix 
  tauhat_pl_mt <- matrix(rep(0,n_test*iter), ncol = iter)
  tauhat_ct_mt <- matrix(rep(0,n_test*iter), ncol = iter)
  tauhat_cf_mt <- matrix(rep(0,n_test*iter), ncol = iter)
  tauhat_xl_mt <- matrix(rep(0,n_test*iter), ncol = iter)
  tauhat_mom_mt <- matrix(rep(0,n_test*iter), ncol = iter)
  true_effect_mt <- matrix(rep(0,n_test*iter), ncol = iter)
  set.seed(2010)
  for(i in 1:iter){
    print(paste("ITER:", i))
    res_d <- simulation_3(d)
    tauhat_pl_mt[,i] <- res_d$tauhat_pl
    tauhat_ct_mt[,i] <- res_d$tauhat_ct
    tauhat_cf_mt[,i] <- res_d$tauhat_cf
    tauhat_xl_mt[,i] <- res_d$tauhat_xl
    tauhat_mom_mt[,i] <- res_d$tauhat_mom
    true_effect_mt[,i] <- res_d$true_effect
  }
  
  # Define comparison metrics for ML estimators
  MSE <- function(tauhat_mt){
    df <- (tauhat_mt - true_effect_mt)^2
    MSE <- apply(df,1,mean)
    return(MSE)
  }
  
  Bias <- function(tauhat_mt){
    Bias <- abs(apply(tauhat_mt-true_effect_mt,1,mean))
    return(Bias)
  }
  
  Coverage <- function(tauhat_mt){
    sehat <- apply(tauhat_mt,1,sd)
    df <- abs(tauhat_mt-true_effect_mt) <= 1.96 * sehat
    Coverage <- apply(df,1,mean)
    return(Coverage)
  }
  
  # Calculate measures
  list_of_ests <- list(tauhat_pl_mt,tauhat_ct_mt,tauhat_cf_mt,tauhat_xl_mt,tauhat_mom_mt)
  
  res_MSE <- sapply(list_of_ests, MSE)
  res_MSE <- as.data.frame(res_MSE)
  colnames(res_MSE) <- c("MSE_pl","MSE_ct","MSE_cf","MSE_xl","MSE_mom")

  res_Bias <- sapply(list_of_ests, Bias)
  res_Bias <- as.data.frame(res_Bias)
  colnames(res_Bias) <- c("Bias_pl","Bias_ct","Bias_cf","Bias_xl","Bias_mom")

  res_Coverage <- sapply(list_of_ests, Coverage)
  res_Coverage <- as.data.frame(res_Coverage)
  colnames(res_Coverage) <- c("Coverage_pl","Coverage_ct","Coverage_cf","Coverage_xl","Coverage_mom")
  res_Coverage

  Metrics <- bind_cols(res_MSE,res_Bias,res_Coverage)
  Metrics <- apply(Metrics,2,mean)
  
  return(Metrics)
})

# Create tables of results
result_raw

# library(data.table)
# result_raw_df <- transpose(data.frame(result_raw)) 
# colnames(result_raw_df) <- names(result_raw[[1]])
# results_table <-  cbind(d=dvals,result_raw_df)

results_table <-  data.frame(cbind(d=dvals, Reduce(rbind, result_raw)))
row.names(results_table) <- c()

results_table_MSE <- results_table %>% select("d","MSE_pl","MSE_ct","MSE_cf","MSE_xl","MSE_mom")
results_table_Bias <- results_table %>% select("d","Bias_pl","Bias_ct","Bias_cf","Bias_xl","Bias_mom")
results_table_Coverage <- results_table %>% select("d","Coverage_pl","Coverage_ct","Coverage_cf","Coverage_xl","Coverage_mom")


xtab_MSE <-  xtable(results_table_MSE, digits = c(0,0,4,4,4,4,4))
print(xtab_MSE, include.rownames = FALSE)

xtab_Bias <-  xtable(results_table_Bias, digits = c(0,0,4,4,4,4,4))
print(xtab_Bias, include.rownames = FALSE)

xtab_Coverage <-  xtable(results_table_Coverage, digits = c(0,0,4,4,4,4,4))
print(xtab_Coverage, include.rownames = FALSE)


#================================| APPLICATIONS|================================
source("application.R")




