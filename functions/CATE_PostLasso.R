# Estimator 1: POST LASSO
#

CATE_PostLasso <- function(Y,W,X,X_test){
  
  df_train = bind_cols(as.data.frame(Y),as.data.frame(W),as.data.frame(X))
  
  # STEP 1: Select variables using Lasso
  ## Get formula
  interactions <- str_c(colnames(as.data.frame(X)), "*W")
  regressors <- str_c(c(colnames(as.data.frame(X)), interactions), collapse = "+")
  fmla_with_interac <- as.formula(str_c("Y~", regressors))
  ## Get Lasso coefficients
  lasso <- glmnetUtils::cv.glmnet(formula = fmla_with_interac, data = df_train)
  beta_xw <- lasso %>% coef()
  ## Get relevant coefficients
  relevant <-  which(beta_xw != 0) %>%  rownames(beta_xw)[.] %>% .[2:length(.)] # Discard intercept
  
  # STEP 2: Apply OLS to the chosen variables
  ## Get formula
  regressors <- str_c(relevant, collapse = "+")
  fmla_with_relevant <- as.formula(str_c("Y~", regressors))
  ## Get OLS estimators
  postlasso <- lm(fmla_with_relevant, data = df_train)
  ## Predict CATE for post-selection Lasso (on test set)
  pl.pred.1 <- predict(postlasso, newdata = data.frame(X_test, W = rep(1,nrow(X_test))))
  pl.pred.0 <- predict(postlasso, newdata = data.frame(X_test, W = rep(0,nrow(X_test))))
  tauhat_pl <- pl.pred.1 - pl.pred.0
  
  return(tauhat_pl)
}