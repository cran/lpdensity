################################################################################
#' Data-driven Bandwidth Selection for Local Polynomial Density Estimation
#'
#' \code{lpbwdensity} implements the bandwidth selector for local polynomial based density (and derivatives) estimation,
#'   proposed in Cattaneo, Jansson and Ma (2017).
#'
#' Companion command: \code{\link{lpdensity}} for local polynomial density estimation.
#'
#' @param data Numeric vector or one dimensional matrix / data frame, the raw data.
#' @param grid Numeric vector or one dimensional matrix / data frame, the grid on which
#'   density is estimated. When set to default, grid points will be chosen as 0.05-0.95
#'   percentiles of the data, with 0.05 step size.
#' @param bwselect String, the method for data-driven bandwidth selection. This option will be
#'   ignored if \code{bw} is provided. Can be (1) \code{"MSE"} (default, mean squared error-optimal
#'   bandwidth selected for each grid point); or (2) \code{"IMSE"} (integrated MSE-optimal bandwidth,
#'   common for all grid points); (3) \code{"ROT"} (rule-of-thumb bandwidth with Gaussian
#'   reference model); and (4) \code{"IROT"} (integrated rule-of-thumb bandwidth with Gaussian
#'   reference model).
#' @param p Integer, the order of the local-polynomial used to construct point
#'   estimates. Should be larger than 0. (Default is 2.)
#' @param v Integer, the derivative to be estimated. Default is \code{1}, which is the density function.
#'   Should be nonnegative.
#' @param kernel String, the kernel function, should be one of \code{"triangular"}, \code{"uniform"} or
#'   \code{"epanechnikov"}.
#' @param Cweights Numeric vector or one dimensional matrix / data frame, the weights used
#'   for counterfactual distribution construction. Should have the same length as sample size.
#'   This option will be ignored if \code{bwselect} is \code{"ROT"} or \code{"IROT"}.
#' @param Pweights Numeric vector or one dimensional matrix / data frame, the weights used
#'   in sampling. Should have the same length as sample size, and nonnegative.
#'   This option will be ignored if \code{bwselect} is \code{"ROT"} or \code{"IROT"}.
#' @param regularize \code{TRUE} (default) or \code{FALSE}, whether the bandwidth should be
#'   regularized. When set to \code{TRUE}, the bandwidth is chosen such that at least 20 + \code{p} + 1
#'   are available locally.
#'
#' @return
#' \item{BW}{A matrix containing (1) \code{grid} (grid points), (2) \code{bw} (bandwidths), and
#'   (3) \code{nh} (effective/local sample sizes).}
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
#' @seealso \code{\link{lpdensity}}.
#'
#' @examples
#' set.seed(42); X <- rnorm(1000)
#' lpbwdensity(X)
#'
#' @export
lpbwdensity <- function(data, grid=c(), bwselect=c("MSE", "IMSE", "ROT", "IROT"),
                        p=c(), v=c(),
                        kernel=c("triangular", "uniform", "epanechnikov"),
                        Cweights=c(), Pweights=c(), regularize=TRUE) {
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

  # bwselect
  if (length(bwselect) == 0) {
    bwselect <- "MSE"
  } else {
    bwselect <- toupper(bwselect[1])
  }
  if (!bwselect%in%c("MSE", "IMSE", "ROT", "IROT")) stop("Incorrect bandwidth selection method specified.\n")

  # p
  if (length(p) == 0) {
    flag_no_p <- TRUE
    p <- 2
  } else if ((length(p) != 1) | !(p[1]%in%0:20)) {
    stop("Polynomial order p incorrectly specified.\n")
  } else {
    flag_no_p <- FALSE
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

  # regularize
  if (length(regularize) == 0) {
    regularize <- FALSE
  } else if (length(regularize) > 1 | !regularize[1]%in%c(TRUE, FALSE)) {
    stop("Regularization parameter incorrectly specified.\n")
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
  # Bandwidth Estimation
  ################################################################################

  if (bwselect == "MSE") {
    bw <- bw_MSE(data=data, grid=grid, p=p, v=v, kernel=kernel, Cweights=Cweights, Pweights=Pweights, regularize=regularize)
  } else if (bwselect == "IMSE") {
    bw <- bw_IMSE(data=data, grid=grid, p=p, v=v, kernel=kernel, Cweights=Cweights, Pweights=Pweights, regularize=regularize)
  } else if (bwselect == "ROT") {
    bw <- bw_ROT(data=data, grid=grid, p=p, v=v, kernel=kernel, regularize=regularize)
  } else {
    bw <- bw_IROT(data=data, grid=grid, p=p, v=v, kernel=kernel, regularize=regularize)
  }

  BW <- cbind(grid, bw)
  colnames(BW) <- c("grid", "bw")
  rownames(BW) <- 1:ng

  nh <- rep(NA, ng)
  for (i in 1:ng) {
    nh[i] <- sum(abs(data-grid[i]) <= BW[i, "bw"])
  }

  BW <- cbind(grid, bw, nh)
  colnames(BW) <- c("grid", "bw", "nh")
  rownames(BW) <- 1:ng

  Result <- list(BW=BW,
                 opt=list(p=p, v=v, kernel=kernel, n=n, ng=ng,
                          bwselect=bwselect,
                          data_min=min(data), data_max=max(data),
                          grid_min=min(grid), grid_max=max(grid)))

  class(Result) <- "CJMlpbwdensity"



  return (Result)
}

################################################################################
#' Internal function.
#'
#' @param x Class \code{CJMlpbwdensity} objects.
#'
#' @keywords internal
#' @export
print.CJMlpbwdensity <- function(x, ...) {

  cat("Call: lpbwdensity\n\n")

  cat(paste("Sample size: ", x$opt$n,  ". Min. observation: ", round(x$opt$data_min, 3), ", Max. observation: ", round(x$opt$data_max, 3), ".\n",   sep=""))
  cat(paste("Grid points: ", x$opt$ng, ". Min. grid: ",        round(x$opt$grid_min, 3), ", Max. grid: ",        round(x$opt$grid_max, 3), ".\n\n", sep=""))

  cat(paste("Polynomial order for point estimation: ", x$opt$p, ".\n", sep=""))
  cat(paste("Order of derivative estimated: ", x$opt$v, ".\n", sep=""))
  cat(paste("Kernel function: ", x$opt$kernel, ".\n", sep=""))
  cat(paste("Bandwidths estimated by: ", x$opt$bwselect, ".\n", sep=""))
  cat("\n")
  cat("Use summary(...) to show bandwidths.\n")
}

################################################################################
#' Internal function.
#'
#' @param object Class \code{CJMlpbwdensity} objects.
#'
#' @keywords internal
#' @export
summary.CJMlpbwdensity <- function(object, ...) {
  x <- object
  colnames(x$BW) <- c("      grid", "        bw", "        nh")

  cat("Call: lpbwdensity\n\n")

  cat(paste("Sample size: ", x$opt$n,  ". Min. observation: ", round(x$opt$data_min, 3), ", Max. observation: ", round(x$opt$data_max, 3), ".\n",   sep=""))
  cat(paste("Grid points: ", x$opt$ng, ". Min. grid: ",        round(x$opt$grid_min, 3), ", Max. grid: ",        round(x$opt$grid_max, 3), ".\n\n", sep=""))

  cat(paste("Polynomial order for point estimation: ", x$opt$p, ".\n", sep=""))
  cat(paste("Order of derivative estimated: ", x$opt$v, ".\n", sep=""))
  cat(paste("Kernel function: ", x$opt$kernel, ".\n", sep=""))
  cat(paste("Bandwidths estimated by: ", x$opt$bwselect, ".\n", sep=""))
  cat("\n")
  cat("Estimation result:\n")
  print(round(x$BW, 3))
}


