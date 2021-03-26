# Estimator 2: (HONEST) CAUSAL TREE 
# Reference:  Athey and Imbens [2016]
#

CATE_CausalTree <- function(Y,W,X,X_test){
  
  df_train = bind_cols(as.data.frame(Y),as.data.frame(W),as.data.frame(X))
  
  # STEP 1: Split the data set
  ## Split training data further into tree training and estimation samples
  df_split <- modelr::resample_partition(df_train, c(train=0.5, estim=0.5))
  df_tr <- as_tibble(df_split$train)
  df_est <- as_tibble(df_split$estim)
  
  # STEP 2: Fit the tree
  ## Get formula
  fmla_tree <- as.formula(str_c("Y~", paste(colnames(as.data.frame(X)), collapse = "+")))
  ## Set parameter
  split.Rule.temp <- "CT"    # The splitting option
  cv.option.temp <- "CT"     # Cross validation option
  split.Honest.temp <- TRUE
  cv.Honest.temp <- TRUE
  split.alpha.temp <- 0.5
  cv.alpha.temp <- 0.5
  split.Bucket.temp <- FALSE
  minisize.temp <- 50
  ## Train honest causal tree
  honest_tree <- causalTree::honest.causalTree(formula = fmla_tree,
                                               data = df_tr,
                                               treatment = df_tr$W,
                                               est_data = df_est,
                                               est_treatment = df_est$W,
                                               split.Rule = split.Rule.temp,
                                               split.Honest = split.Honest.temp,
                                               split.Bucket = split.Bucket.temp,
                                               cv.option = cv.option.temp,  
                                               cv.Honest = cv.Honest.temp,
                                               split.alpha = split.alpha.temp,
                                               cv.alpha = cv.alpha.temp,
                                               minsize = minisize.temp,
                                               HonestSampleSize=length(df_split$estim$idx))
  
  # STEP 3: Cross-validate
  ## Prune honest tree to avoid over-fitting
  opcpid <- which.min(honest_tree$cptable[,4]) # cp minimize xerror
  opcp <- honest_tree$cptable[opcpid,1]
  honest_tree_prune <- prune(honest_tree, opcp)
  
  # STEP 4.a: Predict point estimates and standard errors (on estimation sample) 
  ## Predict point estimates
  #tauhat_ct_est <- predict(honest_tree_prune, newdata = df_est)
  ## Create a factor column 'leaf' indicating leaf assignment
  #num_leaves_est <- length(unique(tauhat_ct_est))
  #df_est$leaf <- factor(tauhat_ct_est, labels = seq(num_leaves_est))
  ## Run regression
  #ols_ct_est <- lm(as.formula("Y ~ 0 + leaf + W:leaf"), data = df_est)
  #ols_ct_summary_est <- summary(ols_ct_est)
  #te_summary_est <- coef(ols_ct_summary_est)[(num_leaves_est+1):(2*num_leaves_est), c("Estimate","Std. Error")]
  #te_summary_est
  
  
  # STEP 4.b: Predict point estimates and standard error (on test set)
  ## Predict CATE for honest tree on test set
  tauhat_ct <- predict(honest_tree_prune, newdata = X_test)
  ## Create a factor column 'leaf' indicating leaf assignment
  #num_leaves <- length(unique(tauhat_ct))
  #df_test$leaf <- factor(tauhat_ct, labels = seq(num_leaves))
  #df_tauhat_ct <- data.frame(tauhat_ct,factor(tauhat_ct, labels = seq(num_leaves))) 
  #colnames(df_tauhat_ct) = c("tauhat_ct","leaf")
  #df_tauhat_ct$leaf <- as.numeric(df_tauhat_ct$leaf)
  ## Run regression
  #ols_ct <- lm(as.formula("Y ~ 0 + leaf + W:leaf"), data = df_test)
  #ols_ct_summary <- summary(ols_ct)
  #te_summary <- coef(ols_ct_summary)[(num_leaves+1):(2*num_leaves), c("Estimate","Std. Error")]
  #te_summary <- te_summary %>% as.data.frame() %>% mutate(., leaf = seq(num_leaves)) %>% rename(.,"tauhat_ct_n" = "Estimate", "sehat_ct" = "Std. Error")
  #df_tauhat_ct <- left_join(df_tauhat_ct,te_summary, by = c("leaf"))
  ## Compute standard errors
  #sehat_ct <- df_tauhat_ct$sehat_ct
  
  return(tauhat_ct)
}