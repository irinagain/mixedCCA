#' From Kendall's tau to correlation matrix estimator using bridge function
#'
#' @param K Kendall's tau correlation matrix (p by p) of \eqn{X}.
#' @param zratio proportion of zero values of dataset \eqn{X}.
#' @param type type of the dataset \eqn{X}. Choose among "continous", "binary", "trunc".
#' @param tol the desired accuracy (convergence tolerance).
#'
#' @examples
#' library(mixedCCA)
#'
#' set.seed(1)
#'
#' n <- 100 # sample size
#' r <- 0.8 # true correlation
#'
#' # Data generation (X: truncated continuous)
#' Z <- mvrnorm(n, mu = c(0, 1, 0), Sigma = matrix(c(1, r, r, r, 1, r, r, r, 1), nrow = 3))
#' X <- ifelse(Z > 0, Z, 0)
#' zratio <- colMeans(X==0)
#' K <- Kendall_matrix(X)
#' fromKtoR(K, zratio = zratio, type = "trunc")
#'
#' @export
fromKtoR <- function(K, zratio = NULL, type = "trunc", tol = 1e-3) {
  K <- as.matrix(K)
  d1 <- nrow(K)
  de1 <- NULL
  if(type == "continuous"){
    hatR <- sin(pi/2 * K)
  }else{
    bridge <- bridge_select(type1 = type, type2 = type)
    ###################################################################
    # Below code is from Yang Ning
    hatR <- matrix(1, d1, d1)

    if ( d1 <= 1 ){
      i <- j <- 1
      f1 <- function(r)(bridge(r, zratio = zratio[i], zratio2 = zratio[j]) - K[i,j])^2
      op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)
      if(op == 100) {
        hatR[i, j] <- hatR[j, i] <-0
      }else {
        hatR[i, j] <- hatR[j, i] <- unlist(op)
      }
    } else {
      for(i in 1:(d1-1)) {
        for(j in (i+1):d1){
          # Below change to use the bridgeF_mix function that was selected previously, no need to supply the type anymore
          f1 <- function(r)(bridge(r, zratio1 = zratio[i], zratio2 = zratio[j]) - K[i,j])^2
          op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)
          if(op == 100) {
            hatR[i, j] <- hatR[j, i] <-0
          }else {
            hatR[i, j] <- hatR[j, i] <- unlist(op)
          }
        }
      }
    }
  }
  return(hatR)
}


#' From Kendall's tau to correlation matrix estimator using bridge function for mixed type data
#'
#' @param K12 Kendall's tau correlation matrix (p1 by p2) between \eqn{X1} and \eqn{X2}.
#' @param zratio1 proportion of zero values of dataset \eqn{X1}.
#' @param zratio2 proportion of zero values of dataset \eqn{X2}.
#' @param type1 type of the dataset \eqn{X1} among "continous", "binary", "trunc".
#' @param type2 type of the dataset \eqn{X2} among "continous", "binary", "trunc".
#' @param tol the desired accuracy (convergence tolerance).
#'
#' @examples
#'
#' library(mixedCCA)
#' set.seed(1)
#'
#' n <- 100 # sample size
#' r <- 0.8 # true correlation
#'
#' # Data generation (X1: truncated continuous, X2: continuous)
#' Z <- mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, r, r, 1), nrow = 2))
#' X1 <- Z[,1]
#' X1[Z[,1] < 1] <- 0
#' X2 <- Z[,2]
#' zratio1 <- mean(X1==0)
#' bridge <- bridge_select(type1 = "trunc", type2 = "continuous")
#' K12 <- KendallTau(X1, X2)
#' # estimated latent correlation R based on K12.
#' fromKtoR_mixed(K12, zratio1 = zratio1, type1 = "trunc", type2 = "continuous")
#'
#' # X1: truncated continuous, X2: truncated continuous
#' X2[Z[,2] < 1] <- 0
#' zratio2 <- mean(X2==0)
#' bridge <- bridge_select(type1 = "trunc", type2 = "trunc")
#' K12 <- KendallTau(X1, X2)
#' K12
#' # estimated latent correlation R based on K12.
#' fromKtoR_mixed(K12, zratio1 = zratio1, zratio2 = zratio2, type1 = "trunc", type2 = "trunc")
#'
#' # X1: truncated continuous, X2: binary
#' X2[X2 != 0] <- 1
#' zratio2 <- mean(X2==0)
#' bridge <- bridge_select(type1 = "trunc", type2 = "binary")
#' K12 <- KendallTau(X1, X2)
#' K12
#' # estimated latent correlation R based on K12.
#' fromKtoR_mixed(K12, zratio1 = zratio1, zratio2 = zratio2, type1 = "trunc", type2 = "binary")
#'
#' @export
fromKtoR_mixed <- function(K12, zratio1 = NULL, zratio2 = NULL, type1 = "trunc", type2 = "continuous", tol = 1e-3) {

  K12 <- as.matrix(K12)
  d1 <- nrow(K12)
  d2 <- ncol(K12)
  ###################################################################
  de1 <- de2 <- NULL

  if(type1 == "continuous" & type2 == "continuous"){
    hatR <- sin(pi/2 * K12)
  }else{
    bridge <- bridge_select(type1 = type1, type2 = type2)
    ###################################################################
    # Below code is from Yang Ning
    hatR <- matrix(1, d1, d2)

    if ( d1 <= 1 & d2 <= 1 ){
      i <- j <- 1
      f1 <- function(r)(bridge(r, zratio1 = zratio1[i], zratio2 = zratio2[j]) - K12[i,j])^2
      op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)
      if(op == 100) {
        hatR[i, j] <- hatR[j, i] <- 0
      }else {
        hatR[i, j] <- hatR[j, i] <- unlist(op)
      }
    } else {
      for(i in 1:d1) {
        for(j in 1:d2){
          # Below change to use the bridgeF_mix function that was selected previously, no need to supply the type anymore
          f1 <- function(r)(bridge(r, zratio1 = zratio1[i], zratio2 = zratio2[j]) - K12[i,j])^2
          op <- tryCatch(optimize(f1, lower = -0.99, upper = 0.99, tol = tol)[1], error = function(e) 100)
          if(op == 100) {
            hatR[i,j] <- 0
          }else {
            hatR[i,j] <- unlist(op)
          }
        }
      }
    }
  }
  return(hatR)
}

