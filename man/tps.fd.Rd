% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tps.fd.R
\name{tps.fd}
\alias{tps.fd}
\title{Tensor product smoothing for functional responses}
\usage{
tps.fd(Y, x, argvals = NULL, k.f = 25, k.x = 15, norder.x = 4,
  basistype.f = "bspline", q.f = 2, q.x = norder.x - 2,
  method = "P-OLS", cov.method = "mod.chol", adapt = FALSE, sp.f = -1,
  sp.x = -1)
}
\arguments{
\item{Y}{\eqn{n \times L} matrix of functional responses}

\item{x}{\eqn{n}-dimensional predictor vector}

\item{argvals}{function argument values}

\item{k.f}{dimension of function-direction basis}

\item{k.x}{dimension of predictor-direction basis}

\item{norder.x}{order of splines in the predictor direction}

\item{basistype.f}{function-direction basis type for the functional direction}

\item{q.x, q.f}{penalty degree for the predictor- and function-direction bases}

\item{method}{"P-OLS" for penalized ordinary least squares, "P-GLS" for penalized 
generalized least squares}

\item{cov.method}{covariance estimation method for GLS: "mod.chol" for modified Cholesky 
decomposition, or "naive"}

\item{adapt}{logical: should the predictor-dimension smoothing parameter be chosen 
in a spatially adaptive manner?}

\item{sp.f, sp.x}{smoothing parameter, passed to \code{\link[mgcv]{bam}}}
}
\value{
An object of class \code{tps.fd}, which is a list whose elements include
\item{coef}{spline coeffients}
\item{fitted}{fitted values}
\item{pwdf}{pointwise degrees of freedom}
}
\description{
Given a set of functional responses and a scalar predictor, this function implements
ordinary least squares and generalized least squares bivariate smooths, with respect to the 
function argument and the predictor
}
\references{
Reiss, P. T., Huang, L., Chen, H., and Colcombe, S. Varying-smoother models for 
functional responses. Available at \url{https://arxiv.org/abs/1412.0778}

Reiss, P. T., Huang, L., Wu, P.-S., Chen, H., and Colcombe, S. Pointwise influence matrices for
functional-response regression. Available at \url{https://works.bepress.com/phil_reiss/43/}
}

