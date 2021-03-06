% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HSVM_2nd_order.R
\name{backtrack_descent}
\alias{backtrack_descent}
\title{Backtracking for steepest descent}
\usage{
backtrack_descent(fx, x, g, d, t = 1, alpha = 0.5, beta = 0.9)
}
\arguments{
\item{fx}{handle to function that returns objective function values}

\item{x}{current parameter estimate}

\item{g}{Gradient vector}

\item{d}{descent direction vector}

\item{t}{current step-size}

\item{alpha}{the backtracking parameter}

\item{beta}{the decrementing multiplier}
}
\description{
Backtracking for steepest descent
}
