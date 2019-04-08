#' Gradient Descent
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
gradient_descent <- function(y, X, b0, delta, t=NULL, lambda=1e-4, max_iter=1e2, tol=1e-3) {
  if(is.numeric(t))
    gradient_descent_fixed(y=y, X=X, b0=b0, delta=delta,t=t,lambda=lambda, max_iter=max_iter, tol=tol)
  else
    gradient_descent_backtrack(y=y, X=X, b0=b0, delta=delta, lambda=lambda, max_iter=max_iter, tol=tol)
}

