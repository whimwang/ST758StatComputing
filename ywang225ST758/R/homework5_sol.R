#' Huberized Hinge Loss
#'
#' @param t input vector
#' @param delta shape parameter
#' @export
#' @examples
#' n <- 1e3
#' delta <- 2
#' t <- seq(-2, 2, length.out=n)
#' y <- phi(t, delta)
#' plot(t, y, type='l')
phi <- function(t, delta) {
  n <- length(t)
  y <- double(n)
  ix <- which(t <= 1 - delta)
  y[ix] <- 1 - t[ix] - delta / 2
  ix <- which(t > 1 - delta & t <= 1)
  y[ix] <- 0.5*((1 - t[ix])^2) / delta
  return(y)
}

#' Huberized Hinge Loss Derivative
#'
#' @param t input vector
#' @param delta shape parameter
#' @export
#' @examples
#' n <- 1e3
#' delta <- 2
#' t <- seq(-2, 2, length.out=n)
#' y <- dphi(t, delta)
#' plot(t, y, type='l')
dphi <- function(t, delta) {
  n <- length(t)
  y <- double(n)
  ix <- which(t <= 1 - delta)
  y[ix] <- -1
  ix <- which(t > 1 - delta & t <= 1)
  y[ix] <- (t[ix] - 1) / delta
  return(y)
}

#' Objective Function for HSVM
#'
#' @param y binary response
#' @param X design matrix
#' @param beta regression coefficient vector
#' @param delta shape parameter
#' @param lambda regularization parameter
#' @export
fx_hsvm <- function(y, X, beta, delta, lambda=1e-4) {
  return(sum(phi(y*(X%*%beta), delta)) + 0.5*lambda*sum(beta ** 2))
}

#' Gradient for HSVM
#'
#' @param y binary response
#' @param X design matrix
#' @param beta regression coefficient vector
#' @param delta shape parameter
#' @param lambda regularization parameter
#' @export
gradf_hsvm <- function(y, X, beta, delta, lambda=1e-4) {
  return(t(X)%*%(y * dphi(y * (X%*%beta), delta)) + lambda*beta)
}

#' Make reduced design matrix
#'
#' @param y binary response
#' @param X design matrix
#' @param beta regression coefficient vector
#' @param delta shape parameter
Xtilde <- function(y, X, beta, delta) {
  t <- y * c(X%*%beta)
  ix <- which(t > 1 - delta & t <= 1)
  return(X[ix,,drop=FALSE])
}

#' Pseudo-Hessian for HSVM
#'
#' @param y binary response
#' @param X design matrix
#' @param beta regression coefficient vector
#' @param delta shape parameter
#' @param lambda regularization parameter
#' @export
dd_hsvm <- function(y, X, beta, delta, lambda) {
  M <- Xtilde(y, X, beta, delta)
  H <- crossprod(M,M)
  diag(H) <- diag(H) + lambda
  return(H)
}

#' Compute Newton Step (Naive) for HSVM
#'
#' @param X Design matrix
#' @param y Binary response vector
#' @param g Gradient vector
#' @param lambda Regularization parameter
#' @param delta shape parameter
#' @export
#' @examples
#' set.seed(12345)
#' n1 <- n2 <- 20; n <- n1 + n2
#' p <- 30
#' mu1 <- sqrt(30)*matrix(abs(rnorm(p)), ncol=1)
#' X <- matrix(rnorm(n*p), n, p)
#' X[1:n1,] <-X[1:n1,] + matrix(rep(mu1, n1), n1, p, byrow=TRUE)
#' y <- rep(-1,n)
#' y[1:n1] <- 1
#' g <- matrix(rnorm(p), ncol=1)
#' beta <- matrix(rnorm(p), ncol=1)
#' lambda <- 1
#' delta <- 1
#' newton_step_naive_sol(X, y, beta, g, lambda, delta)
newton_step_naive_sol <- function(X, y, beta, g, lambda, delta) {
  return(solve(dd_hsvm(y, X, beta, delta, lambda), g))
}

#' Compute Newton Step (Naive) for HSVM
#'
#' @param X Design matrix
#' @param y Binary response vector
#' @param g Gradient vector
#' @param lambda Regularization parameter
#' @param delta shape parameter
#' @export
#' @examples
#' set.seed(12345)
#' n1 <- n2 <- 1e2; n <- n1 + n2
#' p <- 1e3
#' mu1 <- sqrt(30)*matrix(abs(rnorm(p)), ncol=1)
#' X <- matrix(rnorm(n*p), n, p)
#' X[1:n1,] <-X[1:n1,] + matrix(rep(mu1, n1), n1, p, byrow=TRUE)
#' y <- rep(-1,n)
#' y[1:n1] <- 1
#' g <- matrix(rnorm(p), ncol=1)
#' beta <- matrix(rnorm(p), ncol=1)
#' lambda <- 1
#' delta <- 1
#' system.time({db_smw <- newton_step_smw_sol(X, y, beta, g, lambda, delta)})
#' system.time({db_naive <- newton_step_naive_sol(X, y, beta, g, lambda, delta)})
#' plot(db_naive,  db_smw)
#' abline(0,1)
newton_step_smw_sol <- function(X, y, beta, g, lambda, delta) {
  M <- Xtilde(y, X, beta, delta)
  H <- M%*%t(M)
  diag(H) <- diag(H) + lambda*delta
  return((1/lambda) *(g - t(M)%*%solve(H, M%*%g)))
}

#' Backtracking for steepest descent
#'
#' @param fx handle to function that returns objective function values
#' @param x current parameter estimate
#' @param t current step-size
#' @param df the value of the gradient of objective function evaluated at the current x
#' @param d descent direction vector
#' @param alpha the backtracking parameter
#' @param beta the decrementing multiplier
#' @export
#' @examples
#' set.seed(1245)
#' p <- 10
#' fx <- function(x) {
#'    return(0.5*sum(c(x) ** 2))
#' }
#' x0 <- matrix(rnorm(p), ncol=1)
#' df <- x0
#' d <- matrix(0, p, ncol=1)
#' ix <- 1:floor(p/2)
#' d[ix] <- - df[ix]
#' t <- 10
#' backtrack_descent_sol(fx=fx, x=x0, t=t, df=df, d=d)
backtrack_descent_sol <- function(fx, x, t, df, d, alpha=0.5, beta=0.9) {
  f <- fx(x)
  imax <- -14/log10(beta)
  for (i in 1:imax) {
    x_plus <- x + t*d
    f_plus <- fx(x_plus)
    RHS <- f + alpha*t*sum(d*df)
    if (f_plus < RHS) break
    t <- beta*t
  }
  return(t)
}
