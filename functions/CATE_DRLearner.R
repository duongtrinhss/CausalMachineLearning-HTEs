# Estimator 5: DR - LEARNER
# Reference:  Knaus et al. [2021]
#
# Reference: "https://github.com/MCKnaus/CATEs/tree/master/R"

#' This function creates the nuisance parameters p(x), mu(x), and mu_d(x)
#' via cross-fitting using the \code{\link{grf}} package
#'
#' @param y Vector of outcome values
#' @param w Vector of treament indicators
#' @param x Matrix of covariates (N x p matrix)
#' @param index List indicating indices for cross-fitting (e.g. obtained by \code{createFolds} of \code{\link{caret}} pkg)
#' @param args_p List of arguments passed to estimate propensity score model
#' @param args_y List of arguments passed to estimate outcome model
#' @param args_y1 List of arguments passed to estimate outcome model of treated
#' @param args_y0 List of arguments passed to estimate outcome model of non-treated
#' @import grf
#'
#' @return Returns n x 4 matrix containing the nuisance parameters
#'
#' @export

nuisance_cf_grf <- function(y,w,x,index,
                               args_p = list(),
                               args_y =  list(),
                               args_y0 = list(),
                               args_y1 = list()) {
  
  num.trees <- 2000
  np <-  matrix(NA,length(w),4)
  colnames(np) <- c("p_hat","y_hat","y0_hat","y1_hat")
  
  for(i in 1:length(index)) {
    # P-score
    fit_p <- do.call(regression_forest, c(list(X=x[-index[[i]],,drop=F],
                                               Y=w[-index[[i]]]),
                                          #tune.parameters = "none",
                                          num.trees = num.trees,
                                          args_p))
    np[index[[i]],1] <-  predict(fit_p,x[index[[i]],,drop=F])$prediction
    
    # Outcome
    fit_y <- do.call(regression_forest, c(list(X=x[-index[[i]],,drop=F],
                                               Y=y[-index[[i]]]),
                                          #tune.parameters = "none",
                                          num.trees = num.trees,
                                          args_y))
    np[index[[i]],2] <-  predict(fit_y,x[index[[i]],,drop=F])$prediction
    
    # Outcome of non-treated
    fit_y0 <- do.call(regression_forest, c(list(X=x[-index[[i]],,drop=F][w[-index[[i]]]==0,,drop=F],
                                                Y=y[-index[[i]]][w[-index[[i]]]==0]),
                                           #tune.parameters = "none",
                                           num.trees = num.trees,
                                           args_y0))
    np[index[[i]],3] <-  predict(fit_y0,x[index[[i]],,drop=F])$prediction
    
    # Outcome of treated
    fit_y1 <- do.call(regression_forest, c(list(X=x[-index[[i]],,drop=F][w[-index[[i]]]==1,,drop=F],
                                                Y=y[-index[[i]]][w[-index[[i]]]==1]),
                                           #tune.parameters = "none",
                                           num.trees = num.trees,
                                           args_y1))
    np[index[[i]],4] <-  predict(fit_y1,x[index[[i]],,drop=F])$prediction
  }
  
  return(np)
}



#' Implementation of MOM DR using the \code{\link{grf}} package
#'
#' @param y Vector of outcome values
#' @param w Vector of treament indicators
#' @param x Matrix of covariates (N x p matrix) for estimation
#' @param np Matrix of nuisance parameters obtained by by \code{nuisance_cf_glmnet} or \code{nuisance_cf_grf}
#' @param xnew Matrix of covariates (N x p matrix) for out-of-sample prediction
#' @param args_tau List of arguments passed to estimate IATEs
#' @import grf
#'
#' @return Returns vector containing the out-of-sample IATEs
#'
#' @export

mom_dr_grf = function(y,d,x,np,xnew,args_tau=list()) {
  num.trees <- 2000
  mo = np[,"y1_hat"] - np[,"y0_hat"] + d * (y-np[,"y1_hat"]) / np[,"p_hat"] - (1-d) * (y-np[,"y0_hat"]) / (1-np[,"p_hat"])
  fit_tau = do.call(regression_forest,c(list(X=x,Y=mo),
                                        #tune.parameters = "none",
                                        num.trees = num.trees,
                                        args_tau))
  iate = predict(fit_tau,xnew)$prediction
  return(iate)
}

#' This function implements the 50:50 cross-fitting
#'
#' @param est Vector of outcome values
#' @param w Vector of treament indicators
#' @param x Matrix of covariates (N x p matrix)
#' @param index List indicating indices for cross-fitting (e.g. obtained by \code{createFolds} of \code{\link{caret}} pkg)
#' @param args_p List of arguments passed to estimate propensity score model
#' @param args_y List of arguments passed to estimate outcome model
#' @param args_y1 List of arguments passed to estimate outcome model of treated
#' @param args_y0 List of arguments passed to estimate outcome model of non-treated
#'
#' @return Returns n x 4 matrix containing the nuisance parameters
#' 
#' @export


cf_dml1 <- function(est,y,w,x,np,xnew,index,args_tau=list()) {
  
  iate = matrix(0, nrow(xnew), 1)
  
  for(i in 1:length(index)) {
    iate <- iate + 1/length(index) *
      do.call(mom_dr_grf,list(y[index[[i]]],
                       w[index[[i]]],
                       x[index[[i]],,drop=F],
                       np[index[[i]],],
                       xnew,
                       args_tau=args_tau))
  }
  return(iate)
}

#' This function produces the predicted values of "RF MOM DR" estimator 
#'
#' @param Y Vector of training outcome values
#' @param W Vector of training treament indicators
#' @param X Matrix of training covariates (N x p matrix)
#' @param X_test Matrix of validation covariates (N x p matrix)
#' @import grf glmnet caret
#'
#' @return Returns n x 1 matrix containing the IATEs of the validation sample
#'
#' @export

CATE_DRLearner <- function(Y,W,X,X_test){
  
  # Estimate RF nuisance parameters
  index <-  caret::createFolds(Y, k = 2)
  
  np <-  nuisance_cf_grf(Y,W,X,index)
  
  iate <- do.call(cf_dml1,list(mom_dr_grf,Y,W,X,np,X_test,index)) 
  
  return(iate)
}












  
  
  
  