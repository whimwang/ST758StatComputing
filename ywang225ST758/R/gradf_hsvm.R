#' Gradient for HSVM
#'
#' @param y binary response
#' @param X design matrix
#' @param beta regression coefficient vector
#' @param delta shape parameter
#' @param lambda regularization parameter
#' @export
gradf_hsvm <- function(y, X, beta, delta, lambda=1e-4) {
  #y=c(1,-1);beta=c(2,1);delta=0.2;
  #X=matrix(data=c(1,2,3,4),nrow=2,ncol=2,byrow = TRUE);
  p=ncol(X);
  sum_grad=matrix(rep(0,p),nrow=p,ncol=1);
  for(i in 1:nrow(X))
    { 
      t=y[i]*(X[i,]%*%beta);#y_i*x_i_beta
      #print(paste("t value: ",t));
      if(t<=1-delta)
        sum_grad=sum_grad-y[i]*X[i,];
      if((t>1-delta)&(t<=1))
        sum_grad=sum_grad-(y[i]*X[i,])/delta*(1-t);
      #print(paste("sum value: ",n,sum_grad))
   }
  return(sum_grad+lambda*beta)
}


