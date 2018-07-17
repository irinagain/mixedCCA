#' @title Estimate latent correlation matrix based on Kendall's tau
#'
#' @description Under the Gaussian copula model assumption, the latent correlation matrix is estimated based on the observed data of mixed type (continuous/biary/truncated continuous).
#'
#' @param X A numeric data matrix (n by p)
#' @param type A type of data \code{X} among "continuous", "binary", "trunc".
#' @param rho Shrinkage level to make correlation matrix estimator positive definite. The default is 0.01.
#' @return \code{estimateR} returns a data.frame containing
#' \itemize{
#'       \item{type: }{type of the data matrix \code{X}}
#'       \item{R: }{Estimated latent correlation matrix of \code{X}}
#' }
#' @export
#' @importFrom Matrix nearPD
#' @example man/examples/estimateR_ex.R
estimateR <- function(X, type = "trunc", rho = 0.01){
  X <- as.matrix(X)

  n <- nrow(X)
  p1 <- ncol(X)

  if (!(type %in% c("continuous", "binary","trunc"))){
    stop("Unrecognized type of data. Should be one of continuous, binary or trunc.")
  }

  if (type == "continuous"){
    R1 <- sin(pi/2 * pcaPP::cor.fk(X))
  } else {
    zratio1 <- colMeans(X==0)
    R1 <- fromKtoR(Kendall_matrix(X), zratio = zratio1, type = type)
  }

  if ( min(eigen(R1)$values) < 0 ) {
    R1 <- Matrix::nearPD(R1, corr = TRUE)$mat
  }
  # shrinkage method
  R1 <- (1-rho)*R1 + rho*diag(p1)

  return(list(type = type, R = R1))
}
#'
#'
#'
#' @rdname estimateR
#' @aliases estimateR_mixed
#' @param X1 A numeric data matrix (n by p1).
#' @param X2 A numeric data matrix (n by p2).
#' @param type1 A type of data \code{X1} among "continuous", "binary", "trunc".
#' @param type2 A type of data \code{X2} among "continuous", "binary", "trunc".
#' @inheritParams rho
#'
#' @return \code{estimateR_mixed} returns a data.frame containing
#' \itemize{
#'       \item{type1: }{type of the data matrix \code{X1}}
#'       \item{type2: }{type of the data matrix \code{X2}}
#'       \item{R: }{Estimated latent correlation matrix of whole \code{X} = cbind(\code{X1}, \code{X2})}
#'       \item{R1: }{Estimated latent correlation matrix of \code{X1}}
#'       \item{R2: }{Estimated latent correlation matrix of \code{X2}}
#'       \item{R12: }{Estimated latent correlation matrix of \code{X12}}
#' }
#'
#' @export
#' @importFrom Matrix nearPD
estimateR_mixed <- function(X1, X2, type1 = "trunc", type2 = "continuous", rho = 0.01){
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

  if (type[1] == type[2]) {
    # both datasets are of the same type
    if (type[1] == "continuous"){ # both are continuous
      Rall <- sin(pi/2 * pcaPP::cor.fk(cbind(X1, X2)))
      R1 <- Rall[1:p1, 1:p1]
      R2 <- Rall[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)]
      R12 <- Rall[1:p1, (p1 + 1):(p1 + p2)]
    } else {
      zratio1 <- colMeans(X1==0)
      zratio2 <- colMeans(X2==0)
      R1 <- fromKtoR(Kendall_matrix(X1), zratio = zratio1, type = type[1])
      R2 <- fromKtoR(Kendall_matrix(X2), zratio = zratio2, type = type[2])
      R12 <- fromKtoR_mixed(Kendall_matrix(X1, X2), zratio1 = zratio1, zratio2 = zratio2, type1 = type[1], type2 = type[2])
      Rall <- rbind(cbind(R1, R12),cbind(t(R12), R2))
    }
  } else {
    # datasets are of different type
    if (type[1] == "continuous"){
      zratio2 <- colMeans(X2==0)
      R1 <- sin(pi/2 * pcaPP::cor.fk(X1))
      R2 <- fromKtoR(Kendall_matrix(X2), zratio = zratio2, type = type[2])
      R12 <- fromKtoR_mixed(Kendall_matrix(X1, X2), zratio2 = zratio2, type1 = type[1], type2 = type[2])
      Rall <- rbind(cbind(R1, R12), cbind(t(R12), R2))
    } else if (type[2] == "continuous"){
      zratio1 <- colMeans(X1==0)
      R1 <- fromKtoR(Kendall_matrix(X1), zratio = zratio1, type = type[1])
      R12 <- fromKtoR_mixed(Kendall_matrix(X1, X2), zratio1 = zratio1, type1 = type[1], type2 = type[2])
      R2 <- sin(pi/2 * pcaPP::cor.fk(X2))
      Rall <- rbind(cbind(R1, R12), cbind(t(R12), R2))
    } else {
      zratio1 <- colMeans(X1==0)
      zratio2 <- colMeans(X2==0)
      R1 <- fromKtoR(Kendall_matrix(X1), zratio = zratio1, type = type[1])
      R2 <- fromKtoR(Kendall_matrix(X2), zratio = zratio2, type = type[2])
      R12 <- fromKtoR_mixed(Kendall_matrix(X1, X2), zratio1 = zratio1, zratio2 = zratio2, type1 = type[1], type2 = type[2])
      Rall <- rbind(cbind(R1, R12), cbind(t(R12), R2))
    }
  }
  if ( min(eigen(Rall)$values) < 0 ) {
    Rall <- Matrix::nearPD(Rall, corr = TRUE)$mat
  }
  # shrinkage method
  Rall <- (1-rho)*Rall + rho*diag(p1+p2)
  R1 <- Rall[1:p1, 1:p1]
  R2 <- Rall[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)]
  R12 <- Rall[1:p1, (p1 + 1):(p1 + p2)]

  return(list(type = type, R1 = R1, R2 = R2, R12 = R12, R = Rall))
}
