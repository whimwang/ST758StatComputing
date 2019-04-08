context("Homework 5: Huberized Support Vector Machine Part 2")

test_that("Step 1: newton_step_naive is correct", {

  set.seed(12345)
  n1 <- n2 <- 20; n <- n1 + n2
  p <- 30
  mu1 <- sqrt(30)*matrix(abs(rnorm(p)), ncol=1)
  X <- matrix(rnorm(n*p), n, p)
  X[1:n1,] <-X[1:n1,] + matrix(rep(mu1, n1), n1, p, byrow=TRUE)
  y <- rep(-1,n)
  y[1:n1] <- 1
  g <- matrix(rnorm(p), ncol=1)
  beta <- matrix(rnorm(p), ncol=1)
  lambda <- 1
  delta <- 1
  actual <- newton_step_naive(X=X, y=y, b=beta, g=g, lambda=lambda, delta=delta)
  expected <- newton_step_naive_sol(X, y, beta, g, lambda, delta)
  expect_equal( actual, expected )

})

test_that("Step 2: newton_step_smw is correct", {

  set.seed(12345)
  n1 <- n2 <- 20; n <- n1 + n2
  p <- 30
  mu1 <- sqrt(30)*matrix(abs(rnorm(p)), ncol=1)
  X <- matrix(rnorm(n*p), n, p)
  X[1:n1,] <-X[1:n1,] + matrix(rep(mu1, n1), n1, p, byrow=TRUE)
  y <- rep(-1,n)
  y[1:n1] <- 1
  g <- matrix(rnorm(p), ncol=1)
  beta <- matrix(rnorm(p), ncol=1)
  lambda <- 1
  delta <- 1
  actual <- newton_step_smw(X=X, y=y, b=beta, g=g, lambda=lambda, delta=delta)
  expected <- newton_step_smw_sol(X, y, beta, g, lambda, delta)
  expect_equal( actual, expected )

})

test_that("Step 3: backtrack_descent is correct", {

  set.seed(1245)
  p <- 10
  fx <- function(x) {
    return(0.5*sum(c(x) ** 2))
  }
  x0 <- matrix(rnorm(p), ncol=1)
  df <- x0
  d <- matrix(0, p, ncol=1)
  ix <- 1:floor(p/2)
  d[ix] <- - df[ix]
  t <- 10
  actual <- backtrack_descent(fx=fx, x=x0, t=t, d=d, df=df)
  expected <- backtrack_descent_sol(fx=fx, x=x0, t=t, df=df, d=d)
  expect_equal( actual, expected )
})

