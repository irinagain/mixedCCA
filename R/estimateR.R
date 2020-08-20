#' @title Estimate latent correlation matrix
#'
#' @description Estimation of latent correlation matrix from observed data of (possibly) mixed types (continuous/biary/truncated continuous) based on the latent Gaussian copula model.
#'
#' @aliases estimateR estimateR_mixed
#' @param X A numeric data matrix (n by p), n is the sample size and p is the number of variables.
#' @param type A type of variables in \code{X}, must be one of "continuous", "binary" or "trunc".
#' @param method The calculation method of latent correlation. Either "original" method or "approx". If \code{method = "approx"}, multilinear approximation method is used, which is much faster than the original method. If \code{method = "original"}, optimization of the bridge inverse function is used. The default is "approx".
#' @param use.nearPD A logical value indicating whether to use \link[Matrix]{nearPD} or not when the resulting correlation estimator is not positive definite (have at least one negative eigenvalue).
#' @param nu Shrinkage parameter for correlation matrix, must be between 0 and 1, the default value is 0.01.
#' @param tol Desired accuracy when calculating the solution of bridge function.
#' @param verbose If \code{verbose = FALSE}, printing information whether nearPD is used or not is disabled. The default value is FALSE.
#' @return \code{estimateR} returns
#' \itemize{
#'       \item{type: }{Type of the data matrix \code{X}}
#'       \item{R: }{Estimated p by p latent correlation matrix of \code{X}}
#' }
#' @references
#' Fan J., Liu H., Ning Y. and Zou H. (2017) "High dimensional semiparametric latent graphicalmodel for mixed data" <doi:10.1111/rssb.12168>.
#'
#' Yoon G., Carroll R.J. and Gaynanova I. (2020) "Sparse semiparametric canonical correlation analysis for data of mixed types" <doi:10.1093/biomet/asaa007>.
#'
#' Yoon G., Müller C.L., Gaynanova I. (2020) "Fast computation of latent correlations" <arXiv:2006.13875>.
#'
#' @export
#' @import stats
#' @importFrom Matrix nearPD
#' @example man/examples/estimateR_ex.R
estimateR <- function(X, type = "trunc", method = "approx", use.nearPD = TRUE, nu = 0.01, tol = 1e-3, verbose = FALSE){
  X <- as.matrix(X)

  n <- nrow(X)
  p <- ncol(X)

  if (!(type %in% c("continuous", "binary","trunc"))){
    stop("Unrecognized type of data. Should be one of continuous, binary or trunc.")
  }
  ind.NaN <- NULL # initialization of the variable. The default is assuming there is no NA or NaN.

  if (type == "continuous"){
        K <- pcaPP::cor.fk(X)
        R <- sin(pi/2 * K)
  } else if (type == "trunc"){
        # checking data type
        if(sum(X < 0) > 0) {
          stop("The data contains negative values.")
        }
    # checking proportion of zero values
    zratio <- colMeans(X == 0)
        if(sum(zratio) == 0){
          message("The data does not contain zeros. Consider changing the type to \"continuous\".")
        }
        if (sum(zratio == 1) > 0){
          warning("There are variables in the data that have only zeros.\n")
        }

        K <- Kendall_matrix(X)
        if(sum(is.na(K)) > 0){
          warning("There are NaN values in Kendall's tau matrix.\n")
          ind.NaN <- which(colSums(is.na(K)) == (p-1))
          K <- K[-ind.NaN, -ind.NaN]
          zratio <- zratio[-ind.NaN]
        }

        if(method == "approx"){
          R <- fromKtoR_ml(K, zratio = zratio, type = type, tol = tol)
        } else if(method == "original"){
          R <- fromKtoR(K, zratio = zratio, type = type, tol = tol)
        } else {
          stop("Unrecognized method.")
        }

  } else if (type == "binary"){
        # checking data type
        if(sum(!(X %in% c(0, 1))) > 0) {
          stop("The data is not \"binary\".")
        }
        # checking proportion of zero values
        zratio <- colMeans(X == 0)
        if (sum(zratio == 1) > 0 | sum(zratio == 0) > 0){
          warning("There are binary variables in the data that have only zeros or only ones.\n")
        }

        K <- Kendall_matrix(X)
        if(sum(is.na(K)) > 0){
          warning("There are NaN values in Kendall's tau matrix.\n")
          ind.NaN <- which(colSums(is.na(K)) == (p-1))
          K <- K[-ind.NaN, -ind.NaN]
          zratio <- zratio[-ind.NaN]
        }

        if(method == "approx"){
          R <- fromKtoR_ml(K, zratio = zratio, type = type, tol = tol)
        } else if(method == "original"){
          R <- fromKtoR(K, zratio = zratio, type = type, tol = tol)
        } else {
          stop("Unrecognized method.")
        }
  }

  # nearPD
  if ( use.nearPD == TRUE & min(eigen(R)$values) < 0 ) {
        if( verbose ){
          message(" minimum eigenvalue of correlation estimator is ", min(eigen(R)$values), "\n nearPD is used")
        }
        R <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)
  }
  # shrinkage method
  if(nu < 0 | nu > 1){
    stop("nu must be be between 0 and 1.")
  }

  if (length(ind.NaN) > 0){
    R <- (1 - nu)*R + nu*diag(p - length(ind.NaN))
    R.final <- matrix(NaN, nrow = p, ncol = p)
    R.final[-ind.NaN, -ind.NaN] <- R
  } else if (length(ind.NaN) == 0){
    R.final <- (1 - nu)*R + nu*diag(p)
  }

  ### To keep the correct column names for each matrices
  if(length(colnames(X)) == p){
    colnames(R.final) <- rownames(R.final) <- c(colnames(X))
  }

  return(list(type = type, R = R.final))
}

#'
#'
#'
#' @title Estimate latent correlation matrix
#' @rdname estimateR
#' @aliases estimateR estimateR_mixed
#' @param X1 A numeric data matrix (n by p1).
#' @param X2 A numeric data matrix (n by p2).
#' @param type1 A type of variables in \code{X1}, must be one of "continuous", "binary" or "trunc".
#' @param type2 A type of variables in \code{X2}, must be one of "continuous", "binary" or "trunc".
#' @param method The calculation method of latent correlation. Either "original" method or "approx". If \code{method = "approx"}, multilinear approximation method is used, which is much faster than the original method. If \code{method = "original"}, optimization of the bridge inverse function is used. The default is "approx".
#' @param use.nearPD A logical value indicating whether to use \link[Matrix]{nearPD} or not when the resulting correlation estimator is not positive definite (have at least one negative eigenvalue).
#' @param nu Shrinkage parameter for correlation matrix, must be between 0 and 1, the default value is 0.01.
#' @param tol Desired accuracy when calculating the solution of bridge function.
#' @param verbose If \code{verbose = FALSE}, printing information whether nearPD is used or not is disabled. The default value is FALSE.
#'
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
estimateR_mixed <- function(X1, X2, type1 = "trunc", type2 = "continuous", method = "approx", use.nearPD = TRUE, nu = 0.01, tol = 1e-3, verbose = FALSE){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)

  if (nrow(X1) != nrow(X2)){ # Check of they have the same sample size.
    stop ("X1 and X2 must have the same sample size.")
  }

  ind.NaN <- NULL  # initialization of the variable. The default is assuming there is no NA or NaN.

  n <- length(X1)
  p1 <- ncol(X1); p2 <- ncol(X2)

  if (sum(c(type1, type2) %in% c("continuous", "binary", "trunc")) != 2){
    stop("Unrecognised type of variables. Should be one of continuous, binary or trunc.")
  }

  if (type1 == "trunc"){
        if(sum(X1 < 0) > 0) {
          stop("The data X1 contains negative values.")
        }
    zratio1 <- colMeans(X1 == 0)
        if(sum(zratio1) == 0){
          message("The data X1 does not contain zeros. Consider changing the type to \"continuous\".")
        }
        if (sum(zratio1 == 1) > 0){
          warning("There are truncated variables in the data that have only zeros.\n")
        }
  }
  if (type1 == "binary"){
        if(sum(!(X1 %in% c(0, 1))) > 0) {
          stop("The data X1 is not \"binary\".")
        }
    zratio1 <- colMeans(X1 == 0)
        if (sum(zratio1 == 1) > 0 | sum(zratio1 == 0) > 0){
          warning("There are binary variables in the data that have only zeros or only ones.\n")
        }
  }

  if (type2 == "trunc"){
        if(sum(X2 < 0) > 0) {
          stop("The data X2 contains negative values.")
        }
    zratio2 <- colMeans(X2 == 0)
        if(sum(zratio2) == 0){
          message("The data X2 does not contain zeros. Consider changing the type to \"continuous\".")
        }
        if (sum(zratio2 == 1) > 0){
          warning("There are truncated variables in the data that have only zeros.\n")
        }
  }
  if (type2 == "binary"){
        if(sum(!(X2 %in% c(0, 1))) > 0) {
          stop("The data X2 is not \"binary\".")
        }
    zratio2 <- colMeans(X2 == 0)
        if (sum(zratio2 == 1) > 0 | sum(zratio2 == 0) > 0){
          warning("There are binary variables in the data that have only zeros or only ones.\n")
        }
  }

  if (type1 == type2) {
    ################### both datasets are of the same type. CC, TT or BB case.

    Xcomb <- cbind(X1, X2)
    R.final <- estimateR(Xcomb, type = type1, method = method, use.nearPD = use.nearPD, nu = nu, tol = tol)$R
    R1 <- R.final[1:p1, 1:p1]
    R2 <- R.final[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)]
    R12 <- R.final[1:p1, (p1 + 1):(p1 + p2)]

  } else {
    ################### datasets are of different type
    if (type1 == "continuous"){
      ################### These are CT or CB case.
      K1 <- pcaPP::cor.fk(X1)
      R1 <- sin(pi/2 * K1)

      zratio2 <- colMeans(X2 == 0)
      # check Kendall's tau calculation.
      K2 <- Kendall_matrix(X2)
      K12 <- Kendall_matrix(X1, X2)
      if(sum(is.na(K2)) + sum(is.na(K12)) > 0){
        warning("There are NaN values in Kendall's tau matrix.\n")
        ind.NaN <- which(colSums(is.na(K2)) == (p2-1))
        K12 <- K12[, -ind.NaN]
        K2 <- K2[-ind.NaN, -ind.NaN]
        zratio2 <- zratio2[-ind.NaN]
        ind.NaN <- ind.NaN + p1
      }

      if(method == "approx"){
        R2 <- fromKtoR_ml(K2, zratio = zratio2, type = type2, tol = tol)
        R12 <- fromKtoR_ml_mixed(K12, zratio2 = zratio2, type1 = type1, type2 = type2, tol = tol)
      } else if(method == "original"){
        R2 <- fromKtoR(K2, zratio = zratio2, type = type2, tol = tol)
        R12 <- fromKtoR_mixed(K12, zratio2 = zratio2, type1 = type1, type2 = type2, tol = tol)
      } else {
        stop("Unrecognized method.")
      }

      Rall <- rbind(cbind(R1, R12), cbind(t(R12), R2))

    } else if (type2 == "continuous"){
      ################### These are TC or BC case.
      K2 <- pcaPP::cor.fk(X2)
      R2 <- sin(pi/2 * K2)

      zratio1 <- colMeans(X1 == 0)
      # check Kendall's tau calculation.
      K1 <- Kendall_matrix(X1)
      K12 <- Kendall_matrix(X1, X2)
      if(sum(is.na(K1)) + sum(is.na(K12)) > 0){
        warning("There are NaN values in Kendall's tau matrix.\n")
        ind.NaN <- which(colSums(is.na(K1)) == (p1-1))
        K12 <- K12[-ind.NaN, ]
        K1 <- K1[-ind.NaN, -ind.NaN]
        zratio1 <- zratio1[-ind.NaN]
      }

      if(method == "approx"){
        R1 <- fromKtoR_ml(K1, zratio = zratio1, type = type1, tol = tol)
        R12 <- fromKtoR_ml_mixed(K12, zratio1 = zratio1, type1 = type1, type2 = type2, tol = tol)
      } else if(method == "original"){
        R1 <- fromKtoR(K1, zratio = zratio1, type = type1, tol = tol)
        R12 <- fromKtoR_mixed(K12, zratio1 = zratio1, type1 = type1, type2 = type2, tol = tol)
      } else {
        stop("Unrecognized method.")
      }

      Rall <- rbind(cbind(R1, R12), cbind(t(R12), R2))

    } else {
      ################### These are TB or BT case.

      zratio1 <- colMeans(X1 == 0)
      zratio2 <- colMeans(X2 == 0)

      # check Kendall's tau calculation.
      K1 <- Kendall_matrix(X1)
      K2 <- Kendall_matrix(X2)
      K12 <- Kendall_matrix(X1, X2)
      if(sum(is.na(K1)) + sum(is.na(K2)) + sum(is.na(K12)) > 0){
        warning("There are NaN values in Kendall's tau matrix.\n")
        ind.NaN1 <- which(colSums(is.na(K1)) == (p1-1))
        ind.NaN2 <- which(colSums(is.na(K2)) == (p2-1))
        K12 <- K12[-ind.NaN1, -ind.NaN2]
        K1 <- K1[-ind.NaN1, -ind.NaN1]
        K2 <- K2[-ind.NaN2, -ind.NaN2]
        zratio1 <- zratio1[-ind.NaN1]
        zratio2 <- zratio2[-ind.NaN2]
        ind.NaN <- c(ind.NaN1, ind.NaN2+p1)
      }

      if(method == "approx"){
        R1 <- fromKtoR_ml(K1, zratio = zratio1, type = type1, tol = tol)
        R2 <- fromKtoR_ml(K2, zratio = zratio2, type = type2, tol = tol)
        R12 <- fromKtoR_ml_mixed(K12, zratio1 = zratio1, zratio2 = zratio2, type1 = type1, type2 = type2, tol = tol)
      } else if(method == "original"){
        R1 <- fromKtoR(K1, zratio = zratio1, type = type1, tol = tol)
        R2 <- fromKtoR(K2, zratio = zratio2, type = type2, tol = tol)
        R12 <- fromKtoR_mixed(K12, zratio1 = zratio1, zratio2 = zratio2, type1 = type1, type2 = type2, tol = tol)
      } else {
        stop("Unrecognized method.")
      }

      Rall <- rbind(cbind(R1, R12), cbind(t(R12), R2))
    }

    if ( use.nearPD == TRUE & min(eigen(Rall)$values) < 0 ) {
      if( verbose ){
        message(" minimum eigenvalue of correlation estimator is ", min(eigen(Rall)$values), "\n nearPD is used")
      }
      Rall <- as.matrix(Matrix::nearPD(Rall, corr = TRUE)$mat)
    }
    # shrinkage method
    if(nu < 0 | nu > 1){
      stop("nu must be be between 0 and 1.")
    }

    if (length(ind.NaN) > 0){
      Rall <- (1 - nu)*Rall + nu*diag(p1 + p2 - length(ind.NaN))
      R.final <- matrix(NaN, nrow = (p1 + p2), ncol = (p1 + p2))
      R.final[-ind.NaN, -ind.NaN] <- Rall
    } else if (length(ind.NaN) == 0){
      R.final <- (1 - nu)*Rall + nu*diag(p1 + p2)
    }

    ### To keep the correct column names for each matrices
    if(length(colnames(X1)) == p1 & length(colnames(X2)) == p2){
      colnames(R.final) <- rownames(R.final) <- c(colnames(X1), colnames(X2))
    } else if(length(colnames(X1)) != p1 & length(colnames(X2)) == p2){
      colnames(R.final) <- rownames(R.final) <- rep(NA, p1 + p2)
      colnames(R.final)[(p1 + 1):(p1 + p2)] <- rownames(R.final)[(p1 + 1):(p1 + p2)] <- colnames(X2)
    } else if(length(colnames(X1)) == p1 & length(colnames(X2)) != p2){
      colnames(R.final) <- rownames(R.final) <- rep(NA, p1 + p2)
      colnames(R.final)[1:p1] <- rownames(R.final)[1:p1] <- colnames(X1)
    }

    # For convenience, split the R matrices
    R1 <- R.final[1:p1, 1:p1]
    R2 <- R.final[(p1 + 1):(p1 + p2), (p1 + 1):(p1 + p2)]
    R12 <- R.final[1:p1, (p1 + 1):(p1 + p2)]
  }

  return(list(type = c(type1, type2), R1 = R1, R2 = R2, R12 = R12, R = R.final))
}
