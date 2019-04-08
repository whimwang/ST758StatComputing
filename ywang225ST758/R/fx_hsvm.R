#' Objective Function for HSVM
#' 
#' @param y binary response
#' @param X design matrix
#' @param beta regression coefficient vector
#' @param delta shape parameter
#' @param lambda regularization parameter
#' @export
fx_hsvm <- function(y, X, beta, delta, lambda=1e-4) {
  #y=c(1,-1);beta=c(2,1);delta=0.2;
  #X=matrix(data=c(1,2,3,4),nrow=2,ncol=2,byrow = TRUE);
   sum=0;
   for(i in 1:nrow(X))
    {
       t=y[i]*(t(X[i,])%*%beta)
       #print(paste0("t value: ",t))
       if(t<1-delta|t==1-delta)
         sum=sum+(1-t-delta/2);
       if((t>1-delta)&&(t<1|t==1))
        sum=sum+(1-t)^2/2/delta;
        #print(paste0("sum value: ",sum))
     }
  return(sum+lambda/2*sum(beta^2))
}