#' Compute quasi-Newton Step (Naive) for HSVM
#' 
#' @param X Design matrix
#' @param y Binary response vector
#' @param b Regression vector
#' @param g Gradient vector
#' @param delta HSVM parameter
#' @param lambda Regulaization parameter
#' @export
newton_step_naive <- function(X, y, b, g, delta, lambda=1e-4) {
  n=nrow(X);p=ncol(X)
  
  yxb=(X%*%b)*y
  diag_w=ifelse((yxb<=1)&(yxb>1-delta),1/delta,0)
  W=diag(as.vector(diag_w),n)
  
  H=lambda*diag(p)+t(X)%*%W%*%X
  
  R=chol(H)#upper triangle
  y=forwardsolve(t(R),g)
  return(backsolve(R,y)) #for check: all.equal(as.vector(H%*%backsolve(R,y)),g)
}
#' Compute quasi-Newton Step (Sherman-Morrison-Woodbury) for HSVM
#' 
#' @param X Design matrix
#' @param y Binary response vector
#' @param b Regression vector
#' @param g Gradient vector
#' @param delta HSVM parameter
#' @param lambda Regularization parameter
#' @export
newton_step_smw <- function(X, y, b, g, delta, lambda=1e-4) {
  n=nrow(X);p=ncol(X);
  yxb=(X%*%b)*y
  diag_w=ifelse((yxb<=1)&(yxb>1-delta),1/delta,0)
  W=diag(as.vector(diag_w),n)
  
  inverse=diag(p)/lambda-1/(lambda^2)*t(X)%*%solve(diag(n)+W%*%X%*%t(X)/lambda)%*%W%*%X
  return(inverse%*%g) 
  #check:all.equal((lambda*diag(p)+t(X)%*%W%*%X)%*%inverse%*%g,g)
}
#' Backtracking for steepest descent
#' 
#' @param fx handle to function that returns objective function values
#' @param x current parameter estimate
#' @param g Gradient vector
#' @param d descent direction vector
#' @param t current step-size
#' @param alpha the backtracking parameter
#' @param beta the decrementing multiplier
#' @export
backtrack_descent <- function(fx, x, g, d, t=1, alpha=0.5, beta=0.9) {
  t_k=t;
  while(fx(x-t_k*d)>=fx(x)-alpha*t_k*t(g)%*%d)
    t_k=t_k*beta
  return(t_k)
}
#' Damped Newton's Method for Fitting HSVM
#' 
#' @param y Binary response
#' @param X Design matrix
#' @param b Initial regression coefficient vector
#' @param delta HSVM parameter
#' @param lambda regularization parameter
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @param naive which method to use
#' @export
hsvm_newton<- function(X, y, b, delta, lambda=1e-4, max_iter=1e2, tol=1e-3,naive=1) {
  #X=X1;y=y1;b=matrix(rep(0,p),p,1);delta=0.5; lambda=1e-4; max_iter=1e2;tol=1e-3 for debug
  #for backtracking
  x_0=b;
  x_try=x_0;
  times=1;
  while(times<=max_iter)
  {
    f_x=fx_hsvm(y, X, x_try , delta, lambda)
    g_x=gradf_hsvm(y, X, x_try, delta, lambda)
    if(naive==1) d=newton_step_naive(X,y,x_try,g_x,delta,lambda)
    else d=newton_step_smw(X,y,x_try,g_x,delta,lambda)
    #for debug:
    if(t(g_x)%*%d/2<=tol)
      return(x_try)
    step_size=backtrack_descent(fx=function(x) fx_hsvm(y,X,x,delta,lambda),x_try,g_x,d,t=1,alpha=0.5,beta=0.9)
    x_try=x_try-step_size*d
    times<-times+1;
  }
}