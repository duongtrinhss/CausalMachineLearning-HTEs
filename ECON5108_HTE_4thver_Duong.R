#===============================================================================
# Written on April, 2021
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

# Estimator 5: DR - LEARNER
# Reference:  Knaus et al. [2021]
#
source("functions/CATE_DRLearner.R")

#==============| COMPARISON METRICS FOR 4 CAUSAL ML ESTIMATORS |================

#=============================| SIMULATION DESIGN |=============================

# SIMULATION SETUP 1: Low dimensional data; No confounding; 
# Balanced propensity; Linear mean effect; No treatment effect
#
source("simulations/simulation_1.R")

# SIMULATION 2: Low dimensional data; No confounding; Balanced propensity; 
# Nonlinear mean effect; Linear treatment effect
#
source("simulations/simulation_2.R")

# SIMULATION 3: Low dimensional data; No confounding; Balanced propensity; 
# No mean effect; Nonlinear treatment effect
#
# Reference: Wager & Athey(2017)
#
source("simulations/simulation_3.R")

# SIMULATION 4: Low dimensional data; No confounding; Unbalanced propensity propensity; 
# Linear mean effect; Linear treatment effect
#
source("simulations/simulation_4.R")

# SIMULATION 5: Low dimensional data; Confounding; 
# Linear, sparse mean effect; Nonlinear treatment effect
#
source("simulations/simulation_5.R")

# SIMULATION 6: High dimensional data; No Confounding; Balanced propensity propensity; 
# No mean effect; Nonlinear treatment effect
#
source("simulations/simulation_6.R")



#====================| COLLECT RESULTS OF MC SIMULATIONS|=======================
# Parameters for Simulation 1
#
# N = 5000
dvals <- c(5,10,20)
iter <- 30
n_test <- 1000
# simulation: 11111.93 sec elapsed/13305.91 sec elapsed
#
# N = 1000
dvals <- c(5,10,20)
iter <- 150
n_test <- 1000
# simulation: 6626.45 sec elapsed/9471.57 sec

# Parameters for Simulation 2
#
# N = 5000
dvals <- c(5,10,20)
iter <- 30
n_test <- 1000
#
# simulation: 11361.1 sec elapsed
#
# N = 1000
dvals <- c(5,10,20)
iter <- 150
n_test <- 1000
# simulation: 6373.45 sec elapsed/9421.3 sec elapsed

# Parameters for Simulation 3
#
# N = 5000
dvals <- c(2,4,6,8)
iter <- 30
n_test <- 1000
#simulation: 9430.24 sec elapsed
#
# N = 1000
dvals <- c(2,4,6,8)
iter <- 150
n_test <- 1000
# simulation: 5731.19 sec elapsed/8488.31 sec elapsed

# Parameters for Simulation 4
#
# N = 5000
dvals <- c(2,5,8)
iter <- 30
n_test <- 1000
#simulation: 8145.87 sec elapsed
#
# N = 1000
dvals <- c(2,5,8)
iter <- 150
n_test <- 1000
# simulation: 42308.45 sec elapsed

# Parameters for Simulation 5
#
# N = 5000
dvals <- c(2,4,6)
iter <- 30
n_test <- 1000
# simulation: 
#
# N = 1000
dvals <- c(2,4,6)
iter <- 150
n_test <- 1000
# simulation: 5786.2 sec elapsed

# Parameters for Simulation 6
#
dvals <- c(400)
iter <- 150
n_test <- 200
#
dvals <- c(300)
iter <- 30
n_test <- 300
#
dvals <- c(200)
iter <- 30
n_test <- 400
#


tic("simulation")
# Run simulations in turn by replacing line 105
result_raw <- lapply(dvals, function(d){
  print(paste("NOW RUNNING: d = ", d))
  # Get CATEs matrix 
  tauhat_pl_mt <- matrix(rep(0,n_test*iter), ncol = iter)
  tauhat_ct_mt <- matrix(rep(0,n_test*iter), ncol = iter)
  tauhat_cf_mt <- matrix(rep(0,n_test*iter), ncol = iter)
  tauhat_xl_mt <- matrix(rep(0,n_test*iter), ncol = iter)
  tauhat_dr_mt <- matrix(rep(0,n_test*iter), ncol = iter)
  true_effect_mt <- matrix(rep(0,n_test*iter), ncol = iter)
  set.seed(2010)
  for(i in 1:iter){
    print(paste("ITER:", i))
    res_d <- simulation_6(d)
    tauhat_pl_mt[,i] <- res_d$tauhat_pl
    tauhat_ct_mt[,i] <- res_d$tauhat_ct
    tauhat_cf_mt[,i] <- res_d$tauhat_cf
    tauhat_xl_mt[,i] <- res_d$tauhat_xl
    tauhat_dr_mt[,i] <- res_d$tauhat_dr
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
  
  SD <- function(tauhat_mt){
    sehat <- apply(tauhat_mt,1,sd)
    return(sehat)
  }
  
  Coverage <- function(tauhat_mt){
    sehat <- apply(tauhat_mt,1,sd)
    df <- abs(tauhat_mt-true_effect_mt) <= 1.96 * sehat
    Coverage <- apply(df,1,mean)
    return(Coverage)
  }
  
  # Calculate measures
  list_of_ests <- list(tauhat_pl_mt,tauhat_ct_mt,tauhat_cf_mt,tauhat_xl_mt,tauhat_dr_mt)
  
  res_MSE <- sapply(list_of_ests, MSE)
  res_MSE <- as.data.frame(res_MSE)
  colnames(res_MSE) <- c("MSE_pl","MSE_ct","MSE_cf","MSE_xl","MSE_dr")

  res_Bias <- sapply(list_of_ests, Bias)
  res_Bias <- as.data.frame(res_Bias)
  colnames(res_Bias) <- c("Bias_pl","Bias_ct","Bias_cf","Bias_xl","Bias_dr")
  
  res_SD <- sapply(list_of_ests, SD)
  res_SD <- as.data.frame(res_SD)
  colnames(res_SD) <- c("SD_pl","SD_ct","SD_cf","SD_xl","SD_dr")

  res_Coverage <- sapply(list_of_ests, Coverage)
  res_Coverage <- as.data.frame(res_Coverage)
  colnames(res_Coverage) <- c("Coverage_pl","Coverage_ct","Coverage_cf","Coverage_xl","Coverage_dr")
  res_Coverage

  Metrics <- bind_cols(res_MSE,res_Bias,res_SD,res_Coverage)
  Metrics <- apply(Metrics,2,mean)
  
  return(Metrics)
})
toc()
beep(8)



# Create tables of results
result_raw 

library(data.table)
result_raw_df <- transpose(data.frame(result_raw)) 
colnames(result_raw_df) <- names(result_raw[[1]])
results_table <-  cbind(d=dvals,result_raw_df)

results_table <-  data.frame(cbind(d=dvals, Reduce(rbind, result_raw)))
row.names(results_table) <- c()

results_table_MSE <- results_table %>% select("d","MSE_pl","MSE_ct","MSE_cf","MSE_xl","MSE_dr")
results_table_Bias <- results_table %>% select("d","Bias_pl","Bias_ct","Bias_cf","Bias_xl","Bias_dr")
results_table_SD <- results_table %>% select("d","SD_pl","SD_ct","SD_cf","SD_xl","SD_dr")
results_table_Coverage <- results_table %>% select("d","Coverage_pl","Coverage_ct","Coverage_cf","Coverage_xl","Coverage_dr")


xtab_MSE <-  xtable(results_table_MSE, digits = c(0,0,4,4,4,4,4))
print(xtab_MSE, include.rownames = FALSE)

xtab_Bias <-  xtable(results_table_Bias, digits = c(0,0,4,4,4,4,4))
print(xtab_Bias, include.rownames = FALSE)

xtab_SD <-  xtable(results_table_SD, digits = c(0,0,4,4,4,4,4))
print(xtab_SD, include.rownames = FALSE)

xtab_Coverage <-  xtable(results_table_Coverage, digits = c(0,0,4,4,4,4,4))
print(xtab_Coverage, include.rownames = FALSE)


#================================| ADDITIONAL SIMULATION|================================

# Simulation 5b for Causal Forest local centring

# Parameters for Simulation 5b
#
# N = 5000
dvals <- c(2,4,6)
iter <- 30
n_test <- 1000
#simulation: 2208.19 sec elapsed
#
# N = 1000
dvals <- c(2,4,6)
iter <- 150
n_test <- 1000
# simulation: 

tic("simulation")
# Run simulations in turn by replacing line 105
result_raw_add <- lapply(dvals, function(d){
  print(paste("NOW RUNNING: d = ", d))
  # Get CATEs matrix 
  tauhat_cf_lc_mt <- matrix(rep(0,n_test*iter), ncol = iter)
  true_effect_mt <- matrix(rep(0,n_test*iter), ncol = iter)
  set.seed(2010)
  for(i in 1:iter){
    print(paste("ITER:", i))
    res_d <- simulation_5b(d)
    tauhat_cf_lc_mt[,i] <- res_d$tauhat_cf_lc
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
  
  SD <- function(tauhat_mt){
    sehat <- apply(tauhat_mt,1,sd)
    return(sehat)
  }
  
  Coverage <- function(tauhat_mt){
    sehat <- apply(tauhat_mt,1,sd)
    df <- abs(tauhat_mt-true_effect_mt) <= 1.96 * sehat
    Coverage <- apply(df,1,mean)
    return(Coverage)
  }
  
  # Calculate measures
  list_of_ests <- list(tauhat_cf_lc_mt)
  
  res_MSE <- sapply(list_of_ests, MSE)
  res_MSE <- as.data.frame(res_MSE)
  colnames(res_MSE) <- c("MSE_cf_lc")
  
  res_Bias <- sapply(list_of_ests, Bias)
  res_Bias <- as.data.frame(res_Bias)
  colnames(res_Bias) <- c("Bias_cf_lc")
  
  res_SD <- sapply(list_of_ests, SD)
  res_SD <- as.data.frame(res_SD)
  colnames(res_SD) <- c("SD_cf_lc")
  
  res_Coverage <- sapply(list_of_ests, Coverage)
  res_Coverage <- as.data.frame(res_Coverage)
  colnames(res_Coverage) <- c("Coverage_cf_lc")

  Metrics <- bind_cols(res_MSE,res_Bias,res_SD,res_Coverage)
  Metrics <- apply(Metrics,2,mean)
  
  return(Metrics)
})
toc()
beep(8)

result_raw_add


# COMPUTATIONAL TIME

## Include tic()-toc() in each simulation

# print("PostLasso")
# tic()
# tauhat_pl <- CATE_PostLasso(Y,W,X,X_test)
# toc()
# print("CausalTree")
# tic()
# tauhat_ct <- CATE_CausalTree(Y,W,X,X_test)
# toc()
# print("CausalForest")
# tic()
# tauhat_cf <- CATE_CausalForest(Y,W,X,X_test)
# toc()
# print("XLearner")
# tic()
# tauhat_xl <- CATE_XLearner(Y,W,X,X_test)
# toc()
# print("DRLearner")
# tic()
# tauhat_dr <- CATE_DRLearner(Y,W,X,X_test)
# toc()


