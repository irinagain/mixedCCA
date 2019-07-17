#' @title Estimate latent correlation matrix
#'
#'
#' @description Estimation of latent correlation matrix from observed data of (possibly) mixed types (continuous/biary/truncated continuous) based on the latent Gaussian copula model.
#'
#' @aliases estimateR estimateR_mixed
#' @param X A numeric data matrix (n by p), n is the sample size and p is the number of variables.
#' @param type A type of variables in \code{X}, must be one of "continuous", "binary" or "trunc".
#' @param use.nearPD A logical value indicating whether to use \link[Matrix]{nearPD} or not when the resulting correlation estimator is not positive definite (have at least one negative eigenvalue).
#' @param rho Shrinkage parameter for correlation matrix, must be between 0 and 1, the default value is 0.01.
#' @param tol Desired accuracy when calculating the solution of bridge function.
#' @return \code{estimateR} returns
#' \itemize{
#'       \item{type: }{Type of the data matrix \code{X}}
#'       \item{R: }{Estimated p by p latent correlation matrix of \code{X}}
#' }
#' @references
#' Fan J., Liu H., Ning Y. and Zou H. (2017) \href{https://rss.onlinelibrary.wiley.com/doi/abs/10.1111/rssb.12168}{"High dimensional semiparametric latent graphicalmodel for mixed data"}, \emph{J. R. Statist. Soc. B}, 79: 405-421.
#'
#' Yoon G., Carroll R.J. and Gaynanova I. (2018+) \href{http://arxiv.org/abs/1807.05274}{"Sparse semiparametric canonical correlation analysis for data of mixed types"}, \emph{arXiv 1807.05274}.
#' @export
#' @import stats
#' @importFrom Matrix nearPD
#' @example man/examples/estimateR_ex.R
estimateR <- function(X, type = "trunc", use.nearPD = TRUE, rho = 0.01, tol = 1e-3){
  X <- as.matrix(X)

  n <- nrow(X)
  p1 <- ncol(X)

  if (!(type %in% c("continuous", "binary","trunc"))){
    stop("Unrecognized type of data. Should be one of continuous, binary or trunc.")
  }

  if (type == "trunc"){
    if(sum(X<0)>0) {stop("The data contains negative values.")}
    if(sum(colSums(X==0))==0){
      message("The data does not contain zeros. Consider changing the type to \"continuous\".")
    }
  }
  if (type == "binary"){
    if(sum(!(X %in% c(0, 1)))>0) {stop("The data is not \"binary\".")}
  }
  if (type == "continuous"){
    R1 <- sin(pi/2 * pcaPP::cor.fk(X))
  } else {
    zratio1 <- colMeans(X==0)
    R1 <- fromKtoR(Kendall_matrix(X), zratio = zratio1, type = type, tol = tol)
  }

  if ( use.nearPD == TRUE & min(eigen(R1)$values) < 0 ) {
    message(" minimum eigenvalue of correlation estimator is ", min(eigen(R1)$values), "\n nearPD is used")
    R1 <- as.matrix(Matrix::nearPD(R1, corr = TRUE)$mat)
  }
  # shrinkage method
  if(rho<0 | rho>1){ stop("rho must be be between 0 and 1.") }
  R1 <- (1-rho)*R1 + rho*diag(p1)

  return(list(type = type, R = R1))
}
#'
#'
#'
#' @rdname estimateR
#' @title Estimate latent correlation matrix
#'
#' @aliases estimateR estimateR_mixed
#' @param X1 A numeric data matrix (n by p1).
#' @param X2 A numeric data matrix (n by p2).
#' @param type1 A type of variables in \code{X1}, must be one of "continuous", "binary" or "trunc".
#' @param type2 A type of variables in \code{X2}, must be one of "continuous", "binary" or "trunc".
#' @inheritParams use.nearPD
#' @inheritParams rho
#' @inheritParams tol
#'
#' @return \code{estimateR_mixed} returns
#' \itemize{
#'       \item{type1: }{Type of the data matrix \code{X1}}
#'       \item{type2: }{Type of the data matrix \code{X2}}
#'       \item{R: }{Estimated latent correlation matrix of whole \code{X} = (\code{X1}, \code{X2}) (p1+p2 by p1+p2)}
#'       \item{R1: }{Estimated latent correlation matrix of \code{X1} (p1 by p1)}
#'       \item{R2: }{Estimated latent correlation matrix of \code{X2} (p2 by p2)}
#'       \item{R12: }{Estimated latent correlation matrix between \code{X1} and \code{X2} (p1 by p2)}
#' }
#'
#' @export
#' @importFrom Matrix nearPD
estimateR_mixed <- function(X1, X2, type1 = "trunc", type2 = "continuous", use.nearPD = TRUE, rho = 0.01, tol = 1e-3){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)
  if (nrow(X1) != nrow(X2)){ # Check of they have the same length.
    stop ("X1 and X2 must have the same sample size.")
  }
  # if (!(type1 %in% c("continuous", "binary","trunc"))){
  #   stop("Unrecognized type of dataset X1. Should be one of continuous, binary or trunc.")
  # }
  # if (!(type2 %in% c("continuous", "binary","trunc"))){
  #   stop("Unrecognized type of dataset X2. Should be one of continuous, binary or trunc.")
  # }
  n <- length(X1)
  p1 <- ncol(X1); p2 <- ncol(X2)
  type <- rep(NA, 2); type[1] <- type1; type[2] <- type2

  if (sum(type %in% c("continuous", "binary", "trunc")) != 2){
    stop("Unrecognised type of variables. Should be one of continuous, binary or trunc.")
  }

  if (type1 == "trunc"){
    if(sum(X1<0)>0) {stop("The data X1 contains negative values.")}
    if(sum(colSums(X1==0))==0){
      message("The data X1 does not contain zeros. Consider changing the type to \"continuous\".")
    }
  }
  if (type1 == "binary"){
    if(sum(!(X1 %in% c(0, 1)))>0) {stop("The data X1 is not \"binary\".")}
  }

  if (type2 == "trunc"){
    if(sum(X2<0)>0) {stop("The data X2 contains negative values.")}
    if(sum(colSums(X2==0))==0){
      message("The data X2 does not contain zeros. Consider changing the type to \"continuous\".")
    }
  }
  if (type2 == "binary"){
    if(sum(!(X2 %in% c(0, 1)))>0) {stop("The data X2 is not \"binary\".")}
  }

  if (type[1] == type[2]) {
    # both datasets are of the same type
    Xcomb <- cbind(X1, X2)
    Rall <- estimateR(Xcomb, type = type[1], use.nearPD = use.nearPD, rho = rho, tol = tol)$R
    R1 <- Rall[1:p1, 1:p1]
    R2 <- Rall[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)]
    R12 <- Rall[1:p1, (p1 + 1):(p1 + p2)]
  } else {
    # datasets are of different type
    if (type[1] == "continuous"){
      zratio2 <- colMeans(X2==0)
      R1 <- sin(pi/2 * pcaPP::cor.fk(X1))
      R2 <- fromKtoR(Kendall_matrix(X2), zratio = zratio2, type = type[2], tol = tol)
      R12 <- fromKtoR_mixed(Kendall_matrix(X1, X2), zratio2 = zratio2, type1 = type[1], type2 = type[2], tol = tol)
      Rall <- rbind(cbind(R1, R12), cbind(t(R12), R2))
    } else if (type[2] == "continuous"){
      zratio1 <- colMeans(X1==0)
      R1 <- fromKtoR(Kendall_matrix(X1), zratio = zratio1, type = type[1], tol = tol)
      R12 <- fromKtoR_mixed(Kendall_matrix(X1, X2), zratio1 = zratio1, type1 = type[1], type2 = type[2], tol = tol)
      R2 <- sin(pi/2 * pcaPP::cor.fk(X2))
      Rall <- rbind(cbind(R1, R12), cbind(t(R12), R2))
    } else {
      zratio1 <- colMeans(X1==0)
      zratio2 <- colMeans(X2==0)
      R1 <- fromKtoR(Kendall_matrix(X1), zratio = zratio1, type = type[1], tol = tol)
      R2 <- fromKtoR(Kendall_matrix(X2), zratio = zratio2, type = type[2], tol = tol)
      R12 <- fromKtoR_mixed(Kendall_matrix(X1, X2), zratio1 = zratio1, zratio2 = zratio2, type1 = type[1], type2 = type[2], tol = tol)
      Rall <- rbind(cbind(R1, R12), cbind(t(R12), R2))
    }

    if ( use.nearPD == TRUE & min(eigen(Rall)$values) < 0 ) {
      message(" minimum eigenvalue of correlation estimator is ", min(eigen(Rall)$values), "\n nearPD is used")
      Rall <- as.matrix(Matrix::nearPD(Rall, corr = TRUE)$mat)
    }
    # shrinkage method
    Rall <- (1-rho)*Rall + rho*diag(p1+p2)
    R1 <- Rall[1:p1, 1:p1]
    R2 <- Rall[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)]
    R12 <- Rall[1:p1, (p1 + 1):(p1 + p2)]
  }

  return(list(type = type, R1 = R1, R2 = R2, R12 = R12, R = Rall))
}
