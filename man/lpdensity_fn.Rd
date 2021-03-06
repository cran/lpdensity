% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lpdensity_fn.R
\name{lpdensity_fn}
\alias{lpdensity_fn}
\title{Supporting Function for \code{\link{lpdensity}}}
\usage{
lpdensity_fn(
  data,
  grid,
  bw,
  p,
  q,
  v,
  kernel,
  Cweights,
  Pweights,
  massPoints,
  showSE = TRUE
)
}
\arguments{
\item{data}{Numeric vector or one dimensional matrix/data frame, the raw data.}

\item{grid}{Numeric vector or one dimensional matrix/data frame, the grid on which
density is estimated.}

\item{bw}{Numeric vector or one dimensional matrix/data frame, the bandwidth
used for estimation. Should be strictly positive, and have the same length as
\code{grid}.}

\item{p}{Integer, nonnegative, the order of the local-polynomial used to construct point
estimates.}

\item{q}{Integer, nonnegative, the order of the local-polynomial used to construct
confidence interval (a.k.a. the bias correction order).}

\item{v}{Integer, nonnegative, the derivative of distribution function to be estimated. \code{0} for
the distribution function, \code{1} (default) for the density funtion, etc.}

\item{kernel, }{String, the kernel function, should be one of \code{"triangular"},
\code{"uniform"} or \code{"epanechnikov"}.}

\item{Cweights}{Numeric vector or one dimensional matrix/data frame, the weights used
for counterfactual distribution construction. Should have the same length as sample size.}

\item{Pweights}{Numeric vector or one dimensional matrix/data frame, the weights used
in sampling. Should have the same length as sample size, and nonnegative.}

\item{massPoints}{Boolean, whether whether point estimates and standard errors
should be corrected if there are mass points in the data.}

\item{showSE}{\code{TRUE} (default) or \code{FALSE}, whether standard errors should be computed.}
}
\value{
\item{grid}{grid points.}
\item{bw}{bandwidth for each grid point.}
\item{nh}{Effective sample size for each grid point.}
\item{f_p}{Density estimates on the grid with local polynomial of order \code{p},
  with the same length as \code{grid}.}
\item{f_q}{Density estimates on the grid with local polynomial of order \code{q},
  with the same length as \code{grid}. This is reported only if \code{q} is greater than
  0.}
\item{se_p}{Standard errors corresponding to \code{hat_p}.}
\item{se_q}{Standard errors corresponding to \code{hat_q}. This is reported only
  if \code{q} is greater than 0.}
}
\description{
\code{lpdensity_fn} implements the local polynomial density estimator. This
  function is for internal use, and there is no error handling or robustness check.
}
\details{
Recommend: use \code{\link{lpdensity}}.
}
\keyword{internal}
