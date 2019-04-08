#' Softthreshold operator
#'
#' @param x input vector
#' @param lambda threshold
#' @export
#' @examples
#' softthreshold_sol(seq(-10, 10, length.out=10), 1)
softthreshold_sol <- function(x, lambda) {
  return(sign(x)*sapply(abs(x) - lambda, FUN=function(x) {return(max(x, 0))}))
}

#' Proximal mapping of the scaled nuclear norm
#'
#' @param X input matrix
#' @param gamma scale parameter
#' @export
#' @examples
#' set.seed(12345)
#' n <- 10; p <- 20
#' X <- matrix(rnorm(n*p), n, p)
#' Y <- prox_nuc_sol(X, 1)
prox_nuc_sol <- function(X, gamma) {
  S <- svd(X)
  return(S$u %*% diag(softthreshold_sol(S$d, gamma)) %*% t(S$v))
}

#' Project matrix onto observed index set
#'
#' @param X input matrix
#' @param Omega index set
#' @export
#' @examples
#' set.seed(12345)
#' n <- 10; p <- 20
#' X <- matrix(rnorm(n*p), n, p)
#' numel <- n*p
#' Omega <- sample(1:numel, floor(0.5*numel), replace=FALSE)
#' PX <- project_omega_sol(X, Omega)
project_omega_sol <- function(X, Omega) {
  Omegac <- setdiff(1:prod(dim(X)), Omega)
  X[Omegac] <- 0
  return(X)
}

#' Proximal-Gradient Step
#'
#' @param X input matrix
#' @param U current guess
#' @param Omega index set
#' @param lambda regularization parameter
#' @param t step size
#' @examples
#' set.seed(12345)
#' n <- 10; p <- 20
#' X <- matrix(rnorm(n*p), n, p)
#' numel <- n*p
#' Omega <- sample(1:numel, floor(0.5*numel), replace=FALSE)
#' U <- matrix(rnorm(n*p), n, p)
#' U[Omega] <- X[Omega]
#' lambda <- 1
#' t <- 1
#' U1 <- prox_step_sol(X, U, Omega, lambda, t)
prox_step_sol <- function(X, U, Omega, lambda, t) {
  G <- project_omega_sol(U - X, Omega)
  return(prox_nuc_sol(U - t*G, lambda*t))
}

#' Matrix Completion Loss
#'
#' @param X input matrix
#' @param U current guess
#' @param Omega index set
#' @param lambda regularization parameter
#' @param t step size
#' @examples
#' set.seed(12345)
#' n <- 10; p <- 20
#' X <- matrix(rnorm(n*p), n, p)
#' numel <- n*p
#' Omega <- sample(1:numel, floor(0.5*numel), replace=FALSE)
#' U <- matrix(rnorm(n*p), n, p)
#' U[Omega] <- X[Omega]
#' lambda <- 1
#' matrix_completion_loss_sol(X, U, Omega, lambda)
matrix_completion_loss_sol <- function(X, U, Omega, lambda) {
  return(0.5*sum( (X[Omega] - U[Omega]) ** 2) + lambda*sum(svd(U, nu = 0, nv = 0)$d))
}

#' Proximal Gradient Descent
#'
#' @param X input matrix
#' @param U current guess
#' @param Omega index set
#' @param lambda regularization parameter
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @export
#' @examples
#' set.seed(12345)
#' n <- 10; p <- 20
#' X <- matrix(rnorm(n*p), n, p)
#' numel <- n*p
#' Omega <- sample(1:numel, floor(0.5*numel), replace=FALSE)
#' U <- matrix(rnorm(n*p), n, p)
#' U[Omega] <- X[Omega]
#' lambda <- 1
#' pg_sol <- prox_grad_sol(X, U, Omega, lambda)
prox_grad_sol <- function(X, U, Omega, lambda, max_iter=1e2, tol=1e-3) {
  obj <- rel_change_obj <- double(max_iter)
  obj_last <- matrix_completion_loss_sol(X, U, Omega, lambda)
  for (iter in 1:max_iter) {
    U <- prox_step_sol(X, U, Omega, lambda, t=1)
    obj[iter] <- matrix_completion_loss_sol(X, U, Omega, lambda)
    rel_change_obj[iter] <- abs(obj[iter] - obj_last)/(1 + obj_last)
    if (rel_change_obj[iter] < tol) break
    obj_last <- obj[iter]
  }
  return(list(U=U, obj=obj[1:iter], rel_obj=rel_change_obj[1:iter]))
}


#' Soft-Impute Algorithm
#'
#' @param X input matrix
#' @param U current guess
#' @param Omega index set
#' @param lambda regularization parameter
#' @param max_iter maximum number of iterations
#' @param tol convergence tolerance
#' @export
#' @examples
#' set.seed(12345)
#' n <- 10; p <- 20
#' X <- matrix(rnorm(n*p), n, p)
#' numel <- n*p
#' Omega <- sample(1:numel, floor(0.5*numel), replace=FALSE)
#' U <- matrix(rnorm(n*p), n, p)
#' U[Omega] <- X[Omega]
#' lambda <- 1
#' si_sol <- soft_impute_sol(X, U, Omega, lambda)
soft_impute_sol <- function(X, U, Omega, lambda, max_iter=1e2, tol=1e-3) {
  obj <- rel_change_obj <- double(max_iter)
  Omegac <- setdiff(1:prod(dim(X)), Omega)
  M <- X
  obj_last <- matrix_completion_loss_sol(X, U, Omega, lambda)
  for (iter in 1:max_iter) {
    M[Omegac] <- U[Omegac]
    U <- prox_nuc_sol(M, lambda)
    obj[iter] <- matrix_completion_loss_sol(X, U, Omega, lambda)
    rel_change_obj[iter] <- abs(obj[iter] - obj_last)/(1 + obj_last)
    if (rel_change_obj[iter] < tol) break
    obj_last <- obj[iter]
  }
  return(list(U=U, obj=obj[1:iter], rel_obj=rel_change_obj[1:iter]))
}
