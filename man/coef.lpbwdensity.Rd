% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lpbwdensity_methods.R
\name{coef.lpbwdensity}
\alias{coef.lpbwdensity}
\title{Coef Method for Local Polynomial Density Bandwidth Selection}
\usage{
\method{coef}{lpbwdensity}(object, ...)
}
\arguments{
\item{object}{Class "lpbwdensity" object, obtained by calling \code{\link{lpbwdensity}}.}

\item{...}{Other arguments.}
}
\value{
A matrix containing grid points and selected bandwidths.
}
\description{
The coef method for local polynomial density bandwidth selection objects.
}
\examples{
# Generate a random sample
set.seed(42); X <- rnorm(2000)

# Construct bandwidth
coef(lpbwdensity(X))

}
\seealso{
\code{\link{lpbwdensity}} for data-driven bandwidth selection.

Supported methods: \code{\link{coef.lpbwdensity}}, \code{\link{print.lpbwdensity}}, \code{\link{summary.lpbwdensity}}.
}
\author{
Matias D. Cattaneo, Princeton University. \email{cattaneo@princeton.edu}.

Michael Jansson, University of California Berkeley. \email{mjansson@econ.berkeley.edu}.

Xinwei Ma (maintainer), University of California San Diego. \email{x1ma@ucsd.edu}.
}
