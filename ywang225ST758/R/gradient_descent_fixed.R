#' Gradient Descent (Fixed Step-Size)
#'
#' @param y binary response
#' @param X design matrix
#' @param b0 initial parameter estimate
#' @param delta shape parameter
#' @param t step-size
#' @param lambda regularization parameter
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @export
gradient_descent_fixed <- function(y, X, b0, t, delta, lambda=1e-4, max_iter=1e2, tol=1e-3) {
  #y=c(1,-1);beta=c(2,1);delta=0.2;
  #X=matrix(data=c(1,2,3,4),nrow=2,ncol=2,byrow = TRUE);
  print("fixed")
  x_0=b0;times=0;
  while(times<max_iter){
    x_1=x_0-t*gradf_hsvm(y,X,x_0,delta,lambda);
    times=times+1;

      if(abs(fx_hsvm(y, X, beta=x_1, delta, lambda)-fx_hsvm(y, X, beta=x_0, delta, lambda))<tol)# find solution
        return( list(iterate_value=x_1,
                           objective_function_value=fx_hsvm(y, X, beta=x_1, delta, lambda),
                            norm2_gradient_value=sqrt(sum(x_1^2)),
                            relative_change_in_function=abs((fx_hsvm(y, X, beta=x_1, delta, lambda)-fx_hsvm(y, X, beta=x_0, delta, lambda))/fx_hsvm(y,X,beta=x_0,delta,lambda)),
                           relative_change_in_iterate=sqrt(sum((x_1-x_0)^2))/sqrt(sum(x_0^2)) ))  
      else{ x_0=x_1 }# not a solution
  }
    return( list(iterate_value=x_1,
                 objective_function_value=fx_hsvm(y, X, beta=x_1, delta, lambda),
                 norm2_gradient_value=sqrt(sum(x_1^2)),
                 relative_change_in_function=abs((fx_hsvm(y, X, beta=x_1, delta, lambda)-fx_hsvm(y, X, beta=x_0, delta, lambda))/fx_hsvm(y,X,beta=x_0,delta,lambda)),
                 relative_change_in_iterate=sqrt(sum((x_1-x_0)^2))/sqrt(sum(x_0^2)) ))  
}


