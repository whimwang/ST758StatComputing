% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrix_completion_sol.R
\name{prox_step_sol}
\alias{prox_step_sol}
\title{Proximal-Gradient Step}
\usage{
prox_step_sol(X, U, Omega, lambda, t)
}
\arguments{
\item{X}{input matrix}

\item{U}{current guess}

\item{Omega}{index set}

\item{lambda}{regularization parameter}

\item{t}{step size}
}
\description{
Proximal-Gradient Step
}
\examples{
set.seed(12345)
n <- 10; p <- 20
X <- matrix(rnorm(n*p), n, p)
numel <- n*p
Omega <- sample(1:numel, floor(0.5*numel), replace=FALSE)
U <- matrix(rnorm(n*p), n, p)
U[Omega] <- X[Omega]
lambda <- 1
t <- 1
U1 <- prox_step_sol(X, U, Omega, lambda, t)
}
