################################################################################
#' @title lpdensity: Local Polynomial Density Estimation and Inference
#'
#' @description Without imposing stringent distributional assumptions or shape restrictions,
#'   nonparametric density estimation has been popular in economics and other social
#'   sciences for counterfactual analysis, program evaluation, and policy recommendations.
#'   This package implements a novel density estimator based on local polynomial
#'   regression, documented in Cattaneo, Jansson and Ma (2017a): \code{\link{lpdensity}}
#'   to construct local polynomial based density estimator; \code{\link{lpbwdensity}}
#'   to perform data-driven bandwidth selection; and \code{\link{lpdensity.plot}} for
#'   density plot with robust confidence interval.
#'
#'   The companion software article, Cattaneo, Jansson and Ma (2017b),
#'   provides further implementation details and illustrations with simulated data.
#'   For related \code{Stata}
#'   and \code{R} packages useful for nonparametric data analysis  and statistical inference,
#'   visit \url{https://sites.google.com/site/nppackages}.
#'
#' @author
#' Matias D. Cattaneo, University of Michigan. \email{cattaneo@umich.edu}.
#'
#' Michael Jansson, University of California, Berkeley. \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma (maintainer), University of Michigan. \email{xinweima@umich.edu}.
#'
#' @references
#' M. D. Cattaneo, M. Jansson and X. Ma. (2017a). \href{http://www-personal.umich.edu/~cattaneo/papers/Cattaneo-Jansson-Ma_2017_LocPolDensity.pdf}{Simple Local Regression Distribution Estimators}. Working Paper, University of Michigan.
#'
#' M. D. Cattaneo, M. Jansson and X. Ma. (2017b). \href{http://www-personal.umich.edu/~cattaneo/papers/Cattaneo-Jansson-Ma_2017_lpdensity.pdf}{\code{lpdensity}: Local Polynomial Density Estimation and Inference}. Working Paper, University of Michigan.
#'
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom stats qnorm
#' @importFrom stats quantile
#' @importFrom stats D
#' @importFrom stats integrate
#' @importFrom stats optimize
#' @importFrom stats pnorm
#' @importFrom stats dnorm
#' @importFrom stats sd
#' @import ggplot2
#'
#' @aliases lpdensity-package
"_PACKAGE"
