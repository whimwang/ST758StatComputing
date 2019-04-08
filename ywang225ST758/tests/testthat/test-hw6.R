context("Homework 6: Matrix Completion")

test_that("Step 1: softthreshold is correct", {

  n <- 1e3
  x_min <- -1e3
  x_max <- 1e3
  x <- seq(x_min, x_max, length.out=n)
  nLambdas <- 10
  lambdas <- c(0, 10 ** seq(-10, 2*log10(x_max), length.out=nLambdas))

  actual <- expected <- matrix(NA, n, nLambdas)
  for (i in 1:nLambdas) {
    actual[,i] <- softthreshold(x, lambdas[i])
    expected[,i] <- softthreshold_sol(x, lambdas[i])
  }
  expect_equal( actual, expected )

})


test_that("Step 2: prox_nuc is correct", {

  set.seed(1234)
  n <- 20; p <- 30
  X <- matrix(rnorm(n*p), n, p)
  lambda_max <- svd(X, nu=0, nv=0)$d[1]
  nLambdas <- 4
  lambdas <- c(0, 10 ** seq(-10, log10(lambda_max) + 1e-5, nLambdas))

  actual <- expected <- matrix(NA, n*p, nLambdas)
  for (i in 1:nLambdas) {
    actual[,i] <- c(prox_nuc(X, lambdas[i]))
    expected[,i] <- c(prox_nuc_sol(X, lambdas[i]))
  }
  expect_equal( actual, expected )

})

test_that("Step 3: project_omega is correct", {

  set.seed(12345)
  n <- 10; p <- 20
  X <- matrix(rnorm(n*p), n, p)
  numel <- n*p
  Omega <- sample(1:numel, floor(0.5*numel), replace=FALSE)
  actual <- project_omega(X, Omega)
  expected <- project_omega_sol(X, Omega)
  expect_equal(actual, expected)

})

test_that("Step 4: prox_step is correct", {

  set.seed(12345)
  n <- 10; p <- 20
  X <- matrix(rnorm(n*p), n, p)
  numel <- n*p
  Omega <- sample(1:numel, floor(0.5*numel), replace=FALSE)
  U <- matrix(rnorm(n*p), n, p)
  U[Omega] <- X[Omega]
  lambda_max <- svd(U, nu=0, nv=0)$d[1]
  nLambdas <- 4
  lambdas <- c(0, 10 ** seq(-10, log10(lambda_max) + 1e-5, nLambdas))
  nStep <- 3
  steps <- c(0, 1, 2)
  actual <- expected <- matrix(NA, nStep*n*p, nLambdas)
  for (i in 1:nLambdas) {
    for (j in 1:nStep) {
      start <- (j-1)*n*p + 1
      end <- j*n*p
      actual[start:end, i] <- c(prox_step(X, U, Omega, lambdas[i], steps[j]))
      expected[start:end, i] <- c(prox_step_sol(X, U, Omega, lambdas[i], steps[j]))
    }
  }
  expect_equal(actual, expected)
})


test_that("Step 5: matrix_completion_loss is correct", {

  set.seed(12345)
  n <- 10; p <- 20
  X <- matrix(rnorm(n*p), n, p)
  numel <- n*p
  Omega <- sample(1:numel, floor(0.5*numel), replace=FALSE)
  U <- matrix(rnorm(n*p), n, p)
  U[Omega] <- X[Omega]
  lambda_max <- svd(U, nu=0, nv=0)$d[1]
  nLambdas <- 4
  lambdas <- c(0, 10 ** seq(-10, log10(lambda_max) + 1e-5, nLambdas))
  actual <- expected <- matrix(NA, n*p, nLambdas)

  for (i in 1:nLambdas) {
    actual[,i] <- c(matrix_completion_loss(X, U, Omega, lambdas[i]))
    expected[,i] <- c(matrix_completion_loss_sol(X, U, Omega, lambdas[i]))
  }
  expect_equal(actual, expected)

})

