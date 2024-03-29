#' @title Estimate latent correlation matrix
#'
#' @description Estimation of latent correlation matrix from observed data of (possibly) mixed types (continuous/binary/truncated continuous) based on the latent Gaussian copula model.
#'
#' @aliases estimateR estimateR_mixed
#' @param X A numeric data matrix (n by p), n is the sample size and p is the number of variables.
#' @param type A type of variables in \code{X}, must be one of "continuous", "binary" or "trunc".
#' @param method The calculation method of latent correlation. Either "original" method or "approx". If \code{method = "approx"}, multilinear approximation method is used, which is much faster than the original method (requires \code{chebpol} R package). If \code{method = "original"}, optimization of the bridge inverse function is used. The default is "original".
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
#' Yoon G., Mueller C.L., Gaynanova I. (2020) "Fast computation of latent correlations" <arXiv:2006.13875>.
#'
#' @export
#' @import stats
#' @importFrom Matrix nearPD
#' @example man/examples/estimateR_ex.R
estimateR <- function(X, type = "trunc", method = "original", use.nearPD = TRUE, nu = 0.01, tol = 1e-3, verbose = FALSE){
  X <- as.matrix(X)
  p <- ncol(X)

  # shrinkage method
  if(nu < 0 | nu > 1){
    stop("nu must be be between 0 and 1.")
  }

  if (!(type %in% c("continuous", "binary","trunc"))){
    stop("Unrecognized type of data. Should be one of continuous, binary or trunc.")
  }

  ### checking if there is any variable has no variation at all (updated in version 1.4.5)
  ind_sd0 <- which(apply(X, MARGIN = 2, FUN = function(x){ length(unique(x)) == 1 }))
  if(length(ind_sd0) > 0){
    warning("There are variables in the data that have only zeros or only the same values.")
    X <- X[, -ind_sd0, drop = F] # we exclude those variables and only consider this part for rank-based correlation
  }
  if(length(ind_sd0) == p){
    stop("All variables in the data have no variation at all.")
  }

  if (type == "continuous"){
    if (any(is.na(X))){
      # If there are any missing measurements, use slower function
      K <- cor(X, method = "kendall", use = "pairwise.complete.obs")
    }else{
      K <- pcaPP::cor.fk(X)
    }
    R <- sin(pi/2 * K)
  } else {
    zratio <- colMeans(X == 0)
    if (type == "trunc"){
      # checking data type
        if(sum(X < 0) > 0) {
          stop("The data of truncated type contains negative values.")
        }
      # checking proportion of zero values
        if(sum(zratio) == 0){
          message("The data does not contain zeros. Consider changing the type to \"continuous\".")
        }
        # if (sum(zratio == 1) > 0){
        #   warning("There are variables in the data that have only zeros. Filter those variables before continuing. \n")
        # } # deleted in version 1.4.5
    } else {
      # checking data type
      if(sum(!(X %in% c(0, 1))) > 0) {
        stop("The data is not \"binary\".")
      }
      # if (sum(zratio == 1) > 0 | sum(zratio == 0) > 0){
      #   warning("There are binary variables in the data that have only zeros or only ones. Filter those variables before continuing. \n")
      # } # deleted in version 1.4.5
    }
    K <- Kendall_matrix(X)

    if (method == "approx"){
      R <- fromKtoR_ml(K, zratio = zratio, type = type, tol = tol)
    } else {
      R <- fromKtoR(K, zratio = zratio, type = type, tol = tol)
    }
  }

  ### going back to original size of correlation matrix (updated in version 1.4.5)
  if(length(ind_sd0) > 0){
    Rorgsize <- diag(p)
    Rorgsize[-ind_sd0, -ind_sd0] <- R
    R <- Rorgsize
  }

  # nearPD to make it semi pos-definite
  if (use.nearPD == TRUE){
    if (min(eigen(R)$values) < 0) {
      if(verbose){
        message(" minimum eigenvalue of correlation estimator is ", min(eigen(R)$values), "\n nearPD is used")
      }
      R <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)
    }
  }

  # Shrinkage adjustment by nu
  R.final <- (1 - nu) * R + nu * diag(p)

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
estimateR_mixed <- function(X1, X2, type1 = "trunc", type2 = "continuous", method = "original", use.nearPD = TRUE, nu = 0.01, tol = 1e-3, verbose = FALSE){
  X1 <- as.matrix(X1)
  X2 <- as.matrix(X2)

  if (nrow(X1) != nrow(X2)){ # Check of they have the same sample size.
    stop ("X1 and X2 must have the same sample size.")
  }

  # shrinkage method
  if(nu < 0 | nu > 1){
    stop("nu must be be between 0 and 1.")
  }

  p1 <- ncol(X1); p2 <- ncol(X2)

  if (sum(c(type1, type2) %in% c("continuous", "binary", "trunc")) != 2){
    stop("Unrecognised type of variables. Should be one of continuous, binary or trunc.")
  }

  ### checking if there is any variable has no variation at all (updated in version 1.4.5)
  ind1_sd0 <- which(apply(X1, MARGIN = 2, FUN = function(x){ length(unique(x)) == 1 }))
  if(length(ind1_sd0) > 0){
    warning("There are variables in the data X1 that have only zeros or only the same values.")
    X1 <- X1[, -ind1_sd0, drop = F] # we exclude those variables and only consider this part for rank-based correlation
  }
  if(length(ind1_sd0) == p1){
    stop("All variables in the data X1 have no variation at all.")
  }
  ind2_sd0 <- which(apply(X2, MARGIN = 2, FUN = function(x){ length(unique(x)) == 1 }))
  if(length(ind2_sd0) > 0){
    warning("There are variables in the data X2 that have only zeros or only the same values.")
    X2 <- X2[, -ind2_sd0, drop = F] # we exclude those variables and only consider this part for rank-based correlation
  }
  if(length(ind2_sd0) == p2){
    stop("All variables in the data X2 have no variation at all.")
  }

  zratio1 <- colMeans(X1 == 0)
  if (type1 == "trunc"){
    if(sum(X1 < 0) > 0) {
      stop("The data X1 contains negative values.")
    }
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
        if (sum(zratio1 == 1) > 0 | sum(zratio1 == 0) > 0){
          warning("There are binary variables in the data that have only zeros or only ones.\n")
        }
  }

  zratio2 <- colMeans(X2 == 0)
  if (type2 == "trunc"){
    if(sum(X2 < 0) > 0) {
      stop("The data X2 contains negative values.")
    }
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
        if (sum(zratio2 == 1) > 0 | sum(zratio2 == 0) > 0){
          warning("There are binary variables in the data that have only zeros or only ones.\n")
        }
  }

  if (p1 == 1 & p2 == 1){
    # This is just pairwise correlation
    k12 = KendallTau(X1, X2)
    if(method == "approx"){
      r12 = fromKtoR_ml_mixed(k12, zratio1 = zratio1, zratio2 = zratio2, type1 = type1, type2 = type2, tol = tol)
    }else{
      r12 = fromKtoR_mixed(k12, zratio1 = zratio1, zratio2 = zratio2, type1 = type1, type2 = type2, tol = tol)
    }
    return(list(type = c(type1, type2), R1 = 1, R2 = 1, R12 = r12, R = matrix(c(1, r12, r12, 1), 2, 2)))
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
    # Start with 1st dataset
    if (p1 == 1){
      R1 <- 1
    }else if (type1 == "continuous"){
      if (any(is.na(X1))){
        K1 <- cor(X1, method = "kendall", use = "pairwise.complete.obs")
      }else{
        K1 <- pcaPP::cor.fk(X1)
      }
      R1 <- sin(pi/2 * K1)
    }else{
      K1 <- Kendall_matrix(X1)
      if (method == "approx"){
        R1 <- fromKtoR_ml(K1, zratio = zratio1, type = type1, tol = tol)
      } else {
        R1 <- fromKtoR(K1, zratio = zratio1, type = type1, tol = tol)
      }
    }
    # Continue with 2nd dataset
    if (p2 == 1){
      R2 <- 1
    }else if (type2 == "continuous"){
      if (any(is.na(X2))){
        K2 <- cor(X2, method = "kendall", use = "pairwise.complete.obs")
      }else{
        K2 <- pcaPP::cor.fk(X2)
      }
      R2 <- sin(pi/2 * K2)
    }else{
      K2 <- Kendall_matrix(X2)
      if (method == "approx"){
        R2 <- fromKtoR_ml(K2, zratio = zratio2, type = type2, tol = tol)
      } else if(method == "original"){
        R2 <- fromKtoR(K2, zratio = zratio2, type = type2, tol = tol)
      }
    }
    # Do cross-product
    K12 <- Kendall_matrix(X1, X2)

    if (method == "approx"){
      R12 <- fromKtoR_ml_mixed(K12, zratio1 = zratio1, zratio2 = zratio2, type1 = type1, type2 = type2, tol = tol)
    } else {
      R12 <- fromKtoR_mixed(K12, zratio1 = zratio1, zratio2 = zratio2, type1 = type1, type2 = type2, tol = tol)
    }

    Rall <- rbind(cbind(R1, R12), cbind(t(R12), R2))

    ### going back to original size of correlation matrix (updated in version 1.4.5)
    if(length(ind1_sd0) + length(ind2_sd0) > 0){
      Rall_orgsize <- diag(p1 + p2)
      Rall_orgsize[-c(ind1_sd0, p1+ind2_sd0), -c(ind1_sd0, p1+ind2_sd0)] <- Rall
      Rall <- Rall_orgsize
    }

    if (use.nearPD == TRUE){
      if(min(eigen(Rall)$values) < 0) {
        if(verbose) {
          message(" minimum eigenvalue of correlation estimator is ", min(eigen(Rall)$values), "\n nearPD is used")
        }
        Rall <- as.matrix(Matrix::nearPD(Rall, corr = TRUE)$mat)
      }
    }

    # Shrinkage step based on nu
    R.final <- (1 - nu) * Rall + nu * diag(p1 + p2)

    ### To keep the column names in R according to column names that are originally supplied in each matrix
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
