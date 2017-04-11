################################################################################
#' Local Polynomial Density Estimation and Inference
#'
#' \code{lpdensity} implements the local polynomial regression based density (and derivatives)
#'   estimator proposed in Cattaneo, Jansson and Ma (2017). This command can also be
#'   used to obtain smoothed estimates for cumulative distribution functions.
#'
#' Companion command: \code{\link{lpbwdensity}} for data-driven bandwidth selection,
#'   and \code{\link{lpdensity.plot}} for density plot with robust confidence interval.
#'
#' @param data Numeric vector or one dimensional matrix / data frame, the raw data.
#' @param grid Numeric vector or one dimensional matrix / data frame, the grid on which
#'   density is estimated. When set to default, grid points will be chosen as 0.05-0.95
#'   percentiles of the data, with 0.05 step size.
#' @param bw Numeric vector or one dimensional matrix / data frame, the bandwidth
#'   used for estimation. Can be (1) a positive scalar (common
#'   bandwidth for all grid points); or (2) a positive numeric vector specifying bandwidths for
#'   each grid point (should be the same length as \code{grid}).
#' @param p Integer, the order of the local-polynomial used to construct point
#'   estimates. Should be larger than 0. (Default is 2.)
#' @param q Integer, the order of the local-polynomial used to construct pointwise
#'   confidence interval (a.k.a. the bias correction order). Default is \code{p+1}. When specified
#'   \code{0}, no bias correction will be performed. Otherwise it should be larger than 2, and
#'   strictly larger than \code{p}.
#' @param v Integer, the derivative to be estimated. Default is \code{1}, which is the density function.
#'   Should be nonnegative.
#' @param bwselect String, the method for data-driven bandwidth selection. This option will be
#'   ignored if \code{bw} is provided. Can be (1) \code{"MSE"} (default, mean squared error-optimal
#'   bandwidth selected for each grid point); or (2) \code{"IMSE"} (integrated MSE-optimal bandwidth,
#'   common for all grid points); (3) \code{"ROT"} (rule-of-thumb bandwidth with Gaussian
#'   reference model); and (4) \code{"IROT"} (integrated rule-of-thumb bandwidth with Gaussian
#'   reference model).
#' @param kernel String, the kernel function, should be one of \code{"triangular"}, \code{"uniform"} or
#'   \code{"epanechnikov"}.
#' @param Cweights, Numeric vector or one dimensional matrix / data frame, the weights used
#'   for counterfactual distribution construction. Should have the same length as sample size.
#' @param Pweights Numeric vector or one dimensional matrix / data frame, the weights used
#'   in sampling. Should have the same length as sample size and nonnegative.
#' @param scale Numeric, scaling factor for the final estimate. This parameter controls how
#'   estimates are scaled. For example, setting this parameter to 0.5 will scale down both the
#'   point estimates and standard errors by half. By default it is 1. This parameter is used if only
#'   part of the sample is used for estimation, and should not be confused with \code{Cweights}
#'   or \code{Pweights}.
#'
#' @return
#' \item{Estimate}{A matrix containing (1) \code{grid} (grid points), (2) \code{bw} (bandwidths), (3) \code{nh}
#'   (effective/local sample sizes), (4) \code{f_p} (point estimates with p-th order local polynomial),
#'   (5) \code{f_q} (point estimates with q-th order local polynomial, only if option \code{q} is nonzero),
#'   (6) \code{se_p} (standard error corresponding to \code{f_p}), and (7) \code{se_q} (standard error
#'   corresponding to \code{f_q}).}
#' \item{opt}{A list containing options passed to the function.}
#'
#' @references
#' M. D. Cattaneo, M. Jansson and X. Ma. (2017). \href{http://www-personal.umich.edu/~cattaneo/papers/Cattaneo-Jansson-Ma_2017_LocPolDensity.pdf}{Simple Local Regression Distribution Estimators}. Working Paper, University of Michigan.
#'
#' @author
#' Matias D. Cattaneo, University of Michigan. \email{cattaneo@umich.edu}.
#'
#' Michael Jansson, University of California, Berkeley. \email{mjansson@econ.berkeley.edu}.
#'
#' Xinwei Ma (maintainer), University of Michigan. \email{xinweima@umich.edu}.
#'
#' @seealso \code{\link{lpbwdensity}} and \code{\link{lpdensity.plot}}.
#'
#' @examples
#' set.seed(42); X <- rnorm(2000)
#' est1 <- lpdensity(data = X, bwselect = "IMSE")
#' summary(est1)
#'
#' @export
lpdensity <- function(data, grid=c(), bw=c(), p=c(), q=c(), v=c(),
                      bwselect=c("MSE", "IMSE", "ROT", "IROT"),
                      kernel=c("triangular", "uniform", "epanechnikov"),
                      Cweights=c(), Pweights=c(), scale=c()) {
  ################################################################################
  # Input Error Handling
  ################################################################################
  # data
  data <- as.vector(data)
  if (any(is.na(data))) {
    warning("Missing data will be ignored.\n")
    data <- data[!is.na(data)]
  }
  n <- length(data)
  if (!is.numeric(data) | length(data)==0) {
    stop("Data has to be numeric, and cannot be empty.\n")
  }

  # grid
  if (length(grid) == 0) {
    flag_no_grid <- TRUE
    grid <- quantile(data, seq(from=0.05, to=0.95, by=0.05))
    ng <- length(grid)
  } else {
    flag_no_grid <- FALSE
    grid <- as.vector(grid)
    ng <- length(grid)
    if(!is.numeric(grid)) {
      stop("Grid points has to be numeric.\n")
    }
  }

  # bw
  if (length(bw) == 0) {
    flag_no_bw <- TRUE
    if (length(bwselect) == 0) {
      bwselect <- "MSE"
    } else {
      bwselect <- toupper(bwselect[1])
    }
    if (!bwselect%in%c("MSE", "IMSE", "ROT", "IROT")) stop("Incorrect bandwidth selection method specified.\n")
  } else if (length(bw) == 1) {
    if (!is.numeric(bw) | bw <= 0) {
      stop("Bandwidth incorrectly specified.\n")
    } else {
      flag_no_bw <- FALSE
      bw <- rep(bw, ng)
      bwselect <- "MANUAL"
    }
  } else {
    bw <- as.vector(bw)
    if (!is.numeric(bw)) {
      stop("Bandwidth incorrectly specified.\n")
    } else if (length(bw) != ng) {
      stop("Bandwidth has to be the same length as grid.\n")
    } else {
      bwselect <- "MANUAL"
      flag_no_bw <- FALSE
    }
  }

  # p
  if (length(p) == 0) {
    flag_no_p <- TRUE
    p <- 2
  } else if ((length(p) != 1) | !(p[1]%in%0:20)) {
    stop("Polynomial order p incorrectly specified.\n")
  } else {
    flag_no_p <- FALSE
  }

  # q
  if (length(q) == 0) {
    flag_no_q <- TRUE
    q <- p + 1
  } else if ((length(q) > 1) | !(q[1]%in%c(0, 1:20)) | (q[1]<=p & q[1]!=0)) {
    stop("Polynomial order (for bias correction) q incorrectly specified.\n")
  } else {
    flag_no_q <- FALSE
  }

  # v
  if (length(v) == 0) {
    flag_no_v <- TRUE
    v <- 1
  } else if ((length(v) > 1) | !(v[1]%in%c(0:20)) | (v[1]>p)) {
    stop("Derivative order v incorrectly specified.\n")
  } else {
    flag_no_v <- FALSE
  }

  # kernel
  if (length(kernel) == 0) {
    flag_no_kernel <- TRUE
    kernel <- "triangular"
  } else {
    kernel <- tolower(kernel)
    kernel <- kernel[1]
    if (!kernel%in%c("triangular", "uniform", "epanechnikov")) {
      stop("Kernel function incorrectly specified.\n")
    } else {
      flag_no_kernel <- FALSE
    }
  }

  # Cweights
  if (length(Cweights) == 0) {
    flag_no_Cweights <- TRUE
    Cweights <- rep(1, n)
  } else if (!is.numeric(Cweights)) {
    stop("Counterfactual weights incorrectly specified.\n")
  } else if (length(Cweights) != n) {
    stop("Counterfactual weights has to be the same length as sample.\n")
  } else {
    flag_no_Cweights <- FALSE
  }

  # Pweights
  if (length(Pweights) == 0) {
    flag_no_Pweights <- TRUE
    Pweights <- rep(1, n)
  } else if (!is.numeric(Pweights)) {
    stop("Probability weights incorrectly specified.\n")
  } else if (length(Pweights) != n) {
    stop("Probability weights has to be the same length as sample.\n")
  } else if (any(Pweights < 0)) {
    stop("Probability weights has to be nonnegative.\n")
  } else{
    flag_no_Pweights <- FALSE
  }

  # scale
  if (length(scale) == 0) {
    flag_no_scale <- TRUE
    scale <- 1
  } else if (length(scale) > 1 | !is.numeric(scale) | scale < 0) {
    stop("Scale incorrectly specified.\n")
  } else {
    flag_no_scale <- FALSE
  }

  ################################################################################
  # Sample trimming
  ################################################################################
  trim_index <- (Pweights == 0)
  if (all(trim_index)) {
    stop("All weights are zero.\n")
  }
  if (any(trim_index)) {
    data <- data[!trim_index]
    Cweights <- Cweights[!trim_index]
    Pweights <- Pweights[!trim_index]
  }
  if (abs(sum(Cweights * Pweights)) <= .Machine$double.eps * 10) {
    stop("Composited weights (Cweights * Pweights) are numerically zero.\n")
  }
  if (length(data) < 50 + p + 1) stop("Not enough observations.\n")

  ################################################################################
  # Bandwidth selection
  ################################################################################

  if (flag_no_bw & bwselect[1]=="ROT") {
    temp <- bw_ROT(data=data, grid=grid, p=p, v=v, kernel=kernel, regularize=TRUE)
    temp[is.na(temp)] <- 0
    bw <- temp
  }
  if (flag_no_bw & bwselect[1]=="IROT") {
    temp <- bw_IROT(data=data, grid=grid, p=p, v=v, kernel=kernel, regularize=TRUE)
    if (is.na(temp)) temp <- 0
    bw <- rep(temp, ng)
  }
  if (flag_no_bw & bwselect[1]=="MSE") {
    temp <- bw_MSE(data=data, grid=grid, p=p, v=v, kernel=kernel, Cweights=Cweights, Pweights=Pweights, regularize=TRUE)
    temp[is.na(temp)] <- 0
    bw <- temp
  }
  if (flag_no_bw & bwselect[1]=="IMSE") {
    temp <- bw_IMSE(data=data, grid=grid, p=p, v=v, kernel=kernel, Cweights=Cweights, Pweights=Pweights, regularize=TRUE)
    if (is.na(temp)) temp <- 0
    bw <- rep(temp, ng)
  }

  ################################################################################
  # Point Estimation
  ################################################################################
  Estimate <- lpdensity_fn(data=data, grid=grid, bw=bw, p=p, q=q, v=v, kernel=kernel,
                       Cweights=Cweights, Pweights=Pweights)
  Estimate[, c("f_p", "f_q", "se_p", "se_q")] <- Estimate[, c("f_p", "f_q", "se_p", "se_q")] * scale

  row.names(Estimate) <- 1:ng

  Result <- list(Estimate=Estimate,
                 opt=list(
                   p=p, q=q, v=v, kernel=kernel, scale=scale, n=n, ng=ng,
                   bwselect=bwselect,
                   data_min=min(data), data_max=max(data),
                   grid_min=min(grid), grid_max=max(grid)
                 ))
  class(Result) <- "CJMlpdensity"

  ################################################################################
  # Bootstrap
  ################################################################################

  ################################################################################
  # Return
  ################################################################################
  return(Result)
}

################################################################################
#' Internal function.
#'
#' @param x Class \code{CJMlpdensity} objects.
#'
#' @keywords internal
#' @export
print.CJMlpdensity <- function(x, ...) {

  cat("Call: lpdensity\n\n")

  cat(paste("Sample size: ", x$opt$n,  ". Min. observation: ", round(x$opt$data_min, 3), ", Max. observation: ", round(x$opt$data_max, 3), ".\n",   sep=""))
  cat(paste("Grid points: ", x$opt$ng, ". Min. grid: ",        round(x$opt$grid_min, 3), ", Max. grid: ",        round(x$opt$grid_max, 3), ".\n\n", sep=""))

  cat(paste("Polynomial order for point estimation: ", x$opt$p, ".\n", sep=""))
  cat(paste("Order of derivative estimated: ", x$opt$v, ".\n", sep=""))
  if (x$opt$q > 0) {
    cat(paste("Polynomial order for bias-corrected confidence band: ", x$opt$q, ".\n", sep=""))
  } else {
    cat(paste("No bias correction for confidence band.\n", sep=""))
  }
  cat(paste("Kernel function: ", x$opt$kernel, ".\n", sep=""))
  cat(paste("Scaling factor: ",  x$opt$scale,  ".\n", sep=""))
  if (x$opt$bwselect == "MANUAL") {
    cat("Bandwidths provided.\n")
  } else {
    cat(paste("Bandwidths estimated by: ", x$opt$bwselect, ".\n", sep=""))
  }
  cat("\n")
  cat("Use summary(...) to show estimates.\n")
}

################################################################################
#' Internal function.
#'
#' @param object Class \code{CJMlpdensity} objects.
#'
#' @keywords internal
#' @export
summary.CJMlpdensity <- function(object, ...) {
  x <- object
  colnames(x$Estimate) <- c("      grid", "        bw", "        nh",
                            "       f_p", "       f_q", "      se_p", "      se_q")

  cat("Call: lpdensity\n\n")

  cat(paste("Sample size: ", x$opt$n,  ". Min. observation: ", round(x$opt$data_min, 3), ", Max. observation: ", round(x$opt$data_max, 3), ".\n",   sep=""))
  cat(paste("Grid points: ", x$opt$ng, ". Min. grid: ",        round(x$opt$grid_min, 3), ", Max. grid: ",        round(x$opt$grid_max, 3), ".\n\n", sep=""))

  cat(paste("Polynomial order for point estimation: ", x$opt$p, ".\n", sep=""))
  cat(paste("Order of derivative estimated: ", x$opt$v, ".\n", sep=""))
  if (x$opt$q > 0) {
    cat(paste("Polynomial order for bias-corrected confidence band: ", x$opt$q, ".\n", sep=""))
  } else {
    cat(paste("No bias correction for confidence band.\n", sep=""))
  }
  cat(paste("Kernel function: ", x$opt$kernel, ".\n", sep=""))
  cat(paste("Scaling factor: ",  x$opt$scale,  ".\n", sep=""))
  if (x$opt$bwselect == "MANUAL") {
    cat("Bandwidths provided.\n")
  } else {
    cat(paste("Bandwidths estimated by: ", x$opt$bwselect, ".\n", sep=""))
  }
  cat("\n")
  cat("Estimation result:\n")
  print(round(x$Estimate, 3))
}

