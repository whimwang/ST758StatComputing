#' Ridge Regression
#' 
#' \code{ridge_regression} returns the ridge regression coefficient estimates
#' for a sequence of regularization parameter values.
#' 
#' @param y response variables
#' @param X design matrix
#' @param lambda vector of tuning parameters
#' @export
ridge_regression <- function(y, X, lambda) {
  # check conformable
  if(length(y)!=nrow(X)) stop("not conformable");
  # check tuning parameter
  if(length(which(lambda<0))>0) stop("negative tuning parameter exists")
  
  p=ncol(X)
  length_lambda=length(lambda)
  result_beta=matrix(NA,nrow=p,ncol=length_lambda)
  for(i in 1:length_lambda)
   result_beta[,i]=solve(t(X)%*%X+lambda[i]*diag(1,p))%*%t(X)%*%y
  
  return(result_beta) # each col is the beta_est for a lambda
  }