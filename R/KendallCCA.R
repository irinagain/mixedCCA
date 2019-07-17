#' @title Internal data-driven lambda sequence generating function.
#'
#' @description Since this is for sparse CCA, it returns a list of two vectors. Each vector will be used for each data set \eqn{X1} and \eqn{X2}. And \eqn{w1} and \eqn{w2} denote canonical vector for each data set.
#'
#' @param nlamseq The length of lambda sequence
#' @param initlam1 Maximum lambda value for \eqn{X1}
#' @param initlam2 Maximum lambda value for \eqn{X2}
#' @param lam.eps The smallest value for lambda as a fraction of maximum lambda value
#' @param Sigma1 Covariance/correlation matrix of \eqn{X1} (p1 by p1)
#' @param Sigma2 Covariance/correlation matrix of \eqn{X2} (p2 by p2)
#' @param Sigma12 Covariance/correlation matrix betweem \eqn{X1} and \eqn{X2}
#' @param w1init Initial value for canonical vector \eqn{w1}
#' @param w2init Initial value for canonical vector \eqn{w2}
#'
#' @return \code{lambdaseq_generate} returns a list of length 2. Each vector is of the same length \code{nlamseq} and will be used for each data set separately.
#' @export
#'
lambdaseq_generate <- function(nlamseq = 20, initlam1 = NULL, initlam2 = NULL, lam.eps = 1e-02,
                               Sigma1, Sigma2, Sigma12, w1init = NULL, w2init = NULL){
  init.lamseq <- rep(NA, 2)
  p1 <- nrow(Sigma1); p2 <- nrow(Sigma2)
  if (is.null(initlam1)) {
    if (is.null(w2init)) {
      w2init <- rep(1,p2)/sqrt(as.numeric(crossprod(rep(1,p2), Sigma2 %*% rep(1,p2))))
    }
    init.lamseq[1] <- max(abs(Sigma12%*%w2init)) # estimated w1 for different (but close) lamseq could be used
  }

  if (is.null(initlam2)) {
    if (is.null(w1init)) {
      w1init <- rep(1,p1)/sqrt(as.numeric(crossprod(rep(1,p1), Sigma1 %*% rep(1,p1))))
    }
    init.lamseq[2] <- max(abs(t(Sigma12)%*%w1init)) # estimated w2 for different (but close) lamseq could be used
  }

  lamseq <- list("vector")
  # lamseq values in decreasing order.
  lamseq[[1]] <- exp(seq(log(init.lamseq[1]), log(lam.eps*init.lamseq[1]), length.out = nlamseq))
  lamseq[[2]] <- exp(seq(log(init.lamseq[2]), log(lam.eps*init.lamseq[2]), length.out = nlamseq))

  return(lamseq)
}

# wrapper function
#' @title Internal mixedCCA function
#'
#' @description wrapper function for lassobic function in cpp
#'
#' @param n Sample size
#' @param R1 Correlation matrix of dataset \code{X1} (p1 by p1)
#' @param R2 Correlation matrix of dataset \code{X2} (p2 by p2)
#' @param R12 Correlation matrix between the dataset \code{X1} and the dataset \code{X2} (p1 by p2)
#' @param lamseq1 A sequence of lambda values for the datasets \code{X1}. It can be a scalar (a vector of one value). should be the same length with lamseq2.
#' @param lamseq2 A sequence of lambda values for the datasets \code{X2}. It can be a scalar (a vector of one value). should be the same length with lamseq1.
#' @param w1init An initial vector of length p1 for canonical direction \eqn{w1}.
#' @param w2init An initial vector of length p1 for canonical direction \eqn{w2}.
#' @param BICtype Either 1 or 2: For more details for two options, see the reference.
#' @param maxiter The maximum number of iterations allowed.
#' @param tol The desired accuracy (convergence tolerance).
#' @param verbose If \code{verbose = TRUE}, error values in each iteration will be printed. The default value is \code{FALSE}.
#' @param convcheck If \code{convcheck = TRUE}, the convergence error will be printed when the convergence is failed after the algorithm reached \code{maxiter}. The default value is \code{TRUE}.
#'
#' @return \code{find_w12bic} returns a data.frame containing
#' \itemize{
#'       \item{w1: }{estimated canonical direction \eqn{w1}.}
#'       \item{w2: }{estimated canonical direction \eqn{w2}.}
#'       \item{w1trace: }{a matrix, of which column is the estimated canonical direction \eqn{w1} at each iteration. The number of columns is the number of iteration until the convergence.}
#'       \item{w2trace: }{a matrix, of which column is the estimated canonical direction \eqn{w2} at each iteration. The number of columns is the number of iteration until the convergence.}
#'       \item{lam1.iter: }{For each iteration, what lambda value is selected for \eqn{w1} is stored.}
#'       \item{lam2.iter: }{For each iteration, what lambda value is selected for \eqn{w2} is stored.}
#'       \item{obj: }{objective function value without penalty: \eqn{w1^T * R12 * w2}. If lamseq1 and lamseq2 are scalar, then original objective function including penalty part will be used.}
#' }
#' @export
#'
find_w12bic <- function(n, R1, R2, R12, lamseq1, lamseq2, w1init, w2init, BICtype,
                        maxiter = 1000, tol = 0.01, verbose = FALSE, convcheck = TRUE){
  # Same as find_w12 function. SolveLasso part is replaced with lassobic.
  p1 = ncol(R1)
  p2 = ncol(R2)

  error = 1000.0
  iter = 0
  w1 <- w2 <- c()
  wmat1 <- wmat2 <- c()
  lam1.iter <- lam2.iter <- c()
  obj <- c()

  while( iter <= maxiter & error > tol ){
    iter = iter + 1
    error = 0

    ### for w1
    res1 <- lassobic(n = n, R1 = R1, R2 = R2, R12 = R12,
                     w1init = w1init, w2 = w2init, lamseq = lamseq1, BICtype = BICtype,
                     maxiter = maxiter, tol = tol, convcheck = convcheck)
    w1 = res1$finalcoef
    wmat1 <- cbind(wmat1, w1)
    lam1.iter[iter] <- res1$lamseq[res1$bicInd]
    # if w1 is estimated as a vector of only zeros,
    if (sum(abs(w1))==0){
      return(list(w1 = rep(0, p1), w2 = rep(0, p2), lam1.iter = lam1.iter, lam2.iter = lam2.iter, obj = obj))
    }
    w1init = w1

    ### for w2
    res2 <- lassobic(n = n, R1 = R2, R2 = R1, R12 = t(R12),
                     w1init = w2init, w2 = w1init, lamseq = lamseq2, BICtype = BICtype,
                     maxiter = maxiter, tol = tol, convcheck = convcheck)
    w2 = res2$finalcoef
    wmat2 <- cbind(wmat2, w2)
    lam2.iter[iter] <- res2$lamseq[res2$bicInd]
    if (sum(abs(w2))==0){
      return(list(w1 = rep(0, p1), w2 = rep(0, p2), lam1.iter = lam1.iter, lam2.iter = lam2.iter, obj = obj))
    }
    w2init = w2

    if(length(lamseq1) == 1 & length(lamseq2) == 1){
      obj[iter] <- t(w1)%*%R12%*%w2 - lamseq1*sum(abs(w1)) - lamseq2*sum(abs(w2))
    } else {
      obj[iter] <- t(w1)%*%R12%*%w2
    }

    if(iter == 1){
      error <- 1000
    }else if(iter >1){
      error <- abs(obj[iter] - obj[iter-1])/obj[iter-1]
    }

    if (verbose){
      cat("iteration = ", iter, " and selected lambda = ", lam1.iter[iter], "&", lam2.iter[iter], "and objective function = ", obj[iter], "\n")
    }
  }
  if (iter == (maxiter + 1) & convcheck){
    cat("Failed to converge. (findw12 part)\n objective function = ", obj[iter], " and error = ", error, " where tol = ", tol, " with BICtype = ", BICtype, "\n")
  }
  return(list(w1 = w1, w2 = w2, w1trace = wmat1, w2trace = wmat2, lam1.iter = lam1.iter, lam2.iter = lam2.iter, obj = obj))
}




#' @title Sparse CCA for data of mixed types with BIC criterion
#'
#' @description Applies sparse canonical correlation analysis (CCA) for high-dimensional data of mixed types (continuous/biary/truncated continuous). Derived rank-based estimator instead of sample correlation matrix is implemented. There are two types of BIC criteria for variable selection. We found that BIC1 works best for variable selection, whereas BIC2 works best for prediction.
#'
#' @param X1 A numeric data matrix (n by p1).
#' @param X2 A numeric data matrix (n by p2).
#' @param type1 A type of data \code{X1} among "continuous", "binary", "trunc".
#' @param type2 A type of data \code{X2} among "continuous", "binary", "trunc".
#' @param lamseq1 A tuning parameter sequence for \code{X1}. The length should be the same as \code{lamseq2}.
#' @param lamseq2 A tuning parameter sequence for \code{X2}. The length should be the same as \code{lamseq1}.
#' @param initlam1 An initial value to generate a tuning parameter sequence for \code{X1}, which is the maximum value of tuning parameter grid.
#' @param initlam2 An initial value to generate a tuning parameter sequence for \code{X2}, which is the maximum value of tuning parameter grid.
#' @param nlamseq The number of tuning parameter sequence lambda - default is 20.
#' @param lam.eps A ratio of the smallest value for lambda to the maximum value of lambda.
#' @param w1init An initial vector of length p1 for canonical direction \eqn{w1}.
#' @param w2init An initial vector of length p2 for canonical direction \eqn{w2}.
#' @param BICtype Either 1 or 2: For more details for two options, see the reference.
#' @param KendallR An estimated Kendall \eqn{\tau} matrix. The default is NULL, which means that it will be automatically estimated by Kendall's \eqn{\tau} estimator unless the user supplies.
#' @param tol The desired accuracy (convergence tolerance).
#' @param maxiter The maximum number of iterations allowed.
#' @param verbose If \code{verbose = TRUE}, error values in each iteration will be printed. The default value is \code{FALSE}.
#' @param convcheck If \code{convcheck = TRUE}, the convergence error will be printed when the convergence is failed after the algorithm reached \code{maxiter}. The default value is \code{TRUE}.
#'
#' @references
#' Yoon G., Carroll R.J. and Gaynanova I. (2018+) \href{http://arxiv.org/abs/1807.05274}{"Sparse semiparametric canonical correlation analysis for data of mixed types"}.
#' @return \code{mixedCCA} returns a data.frame containing
#' \itemize{
#'       \item{KendallR: }{estimated Kendall's \eqn{\tau} matrix estimator.}
#'       \item{lambda_seq: }{the values of \code{lamseq} used for sparse CCA.}
#'       \item{w1: }{estimated canonical direction \eqn{w1}.}
#'       \item{w2: }{estimated canonical direction \eqn{w2}.}
#'       \item{cancor: }{estimated canonical correlation.}
#'       \item{selected_x1: }{indices of selected variables in \code{X1}.}
#'       \item{selected_x2: }{indices of selected variables in \code{X2}.}
#'       \item{fitresult: }{more details on selected lambda values in each iteration.}
#' }
#'
#' @example man/examples/mixedCCA_ex.R
#' @useDynLib mixedCCA, .registration = TRUE
#' @import stats
#' @importFrom Rcpp evalCpp
#' @export
mixedCCA <- function(X1, X2, type1, type2, lamseq1 = NULL, lamseq2 = NULL, initlam1 = NULL, initlam2 = NULL,
                     nlamseq = 20, lam.eps = 1e-02,
                     w1init = NULL, w2init = NULL, BICtype,
                     KendallR = NULL,
                     tol = 0.01, maxiter = 1000, verbose = FALSE, convcheck = TRUE){
  n <- nrow(X1)
  p1 <- ncol(X1); p2 <- ncol(X2);
  p <- p1 + p2

  ### Compute Kendall tau.
  if(is.null(KendallR)){
    R <- estimateR_mixed(X1, X2, type1 = type1, type2 = type2)$R
  } else {
    R <- KendallR; rm(KendallR)
  }

  R1 <- as.matrix(R[1:p1, 1:p1])
  R2 <- as.matrix(R[(p1+1):p, (p1+1):p])
  R12 <- as.matrix(R[1:p1, (p1+1):p])

  if (is.null(w1init) | is.null(w2init)){
      RCCA <- standardCCA(S1 = R1, S2 = R2, S12 = R12)
      if (is.null(w1init)){
        w1init <- as.matrix(RCCA$w1, ncol=1)
      }
      if (is.null(w2init)){
        w2init <- as.matrix(RCCA$w2, ncol=1)
      }
  }

  ### Create lambda sequences
  lambda_seq <- list()
  if (is.null(lamseq1) & is.null(lamseq2)){
    # Perform standardization of w1 and w2
    w1init <- w1init/sqrt(as.numeric(crossprod(w1init, R1 %*% w1init)))
    w2init <- w2init/sqrt(as.numeric(crossprod(w2init, R2 %*% w2init)))
    lambda_seq <- lambdaseq_generate(nlamseq = nlamseq, initlam1 = initlam1, initlam2 = initlam2, lam.eps = lam.eps,
                                     Sigma1 = R1, Sigma2 = R2, Sigma12 = R12, w1init = w1init, w2init = w2init)
  } else {
    if(length(lamseq1) == length(lamseq2)){
      lambda_seq[[1]] <- lamseq1
      lambda_seq[[2]] <- lamseq2
    } else {
      warning( "The lengths of lambda sequences for two variables are different." );
      stop;
    }
  }

  ### Calculate canonical coefficients using BIC.
  fitresult <- find_w12bic(n = n, R1 = R1, R2 = R2, R12 = R12,
                           lamseq1 = lambda_seq[[1]], lamseq2 = lambda_seq[[2]],
                           w1init = w1init, w2init = w2init, BICtype = BICtype, maxiter = maxiter, tol = tol,
                           verbose = verbose, convcheck = convcheck)

  w1 <- fitresult$w1
  w2 <- fitresult$w2
  ix1 <- which(w1!=0)
  ix2 <- which(w2!=0)

  cancor <- t(w1)%*%R12%*%w2

  return(list(KendallR = R,
              lambda_seq = lambda_seq,
              w1 = w1, w2 = w2,
              cancor = cancor,
              selected_x1 = ix1, selected_x2 = ix2, fitresult = fitresult))
}


# modified to only deal with positive eigenvalues larger than the tolerance.
# inputs are correlation/covariance matrices.
# only returning first pair of canonical covaraites.

#' @title Internal standard CCA function.
#'
#' @description This function is modified from original CCA function to only deal with positive eigenvalues larger than the tolerance. Only returning the first pair of canonical covariates.
#'
#' @param S1 correlation/covariance matrix of dataset \code{X1}.
#' @param S2 correlation/covariance matrix of dataset \code{X2}.
#' @param S12 correlation/covariance matrix between dataset \code{X1} and dataset \code{X2}.
#' @param tol tolerance for eigenvalues. \code{standardCCA} function only deals with positive eigenvalues larger than the tolerance.
#'
#' @return \code{standardCCA} returns a data.frame containing
#' \itemize{
#'       \item{cancor: }{estimated canonical correlation.}
#'       \item{w1: }{estimated canonical direction \eqn{w1}.}
#'       \item{w2: }{estimated canonical direction \eqn{w2}.}
#' }
#'
#' @importFrom irlba irlba
#' @export
#'
standardCCA <- function(S1, S2, S12, tol = 1e-4){
  S1 <- as.matrix(S1)
  S2 <- as.matrix(S2)
  S12 <- matrix(S12, nrow = dim(S1)[1], ncol = dim(S2)[1])

  S1eig <- eigen(S1, symmetric=TRUE)
  subset <- which(S1eig$values > tol)
  S1isqrt <- as.matrix(S1eig$vectors[, subset] %*% diag(1/sqrt(S1eig$values[subset])) %*% t(S1eig$vectors[, subset]))
  S2eig <- eigen(S2, symmetric=TRUE)
  subset <- which(S2eig$values > tol)
  S2isqrt <- as.matrix(S2eig$vectors[, subset] %*% diag(1/sqrt(S2eig$values[subset])) %*% t(S2eig$vectors[, subset]))

  outsvd <- irlba::irlba(S1isqrt %*% S12 %*% S2isqrt, nv = 1) # svd2(S1isqrt %*% S12 %*% S2isqrt, nu = 1, nv = 1)

  # We are considering only one pair of canonical covariates.
  cancor <- max(outsvd$d)

  w1 <- S1isqrt %*% outsvd$u
  w2 <- S2isqrt %*% outsvd$v

  return(list(cancor = cancor, w1 = w1, w2 = w2))
}

# modified from CCA package rcc function
# this estimates are used for initial values for our algorithm. warm starts.
# inputs are rank-based estimator.
# lambda 1 and lambda 2 are scalars.

#' @title Internal RidgeCCA function
#' @description This function is modified from CCA package rcc function. This estimates are used for initial values for our algorithm as warm starts.
#'
#' @param R1 correlation/covariance/rank-based correlation matrix of dataset \code{X1}.
#' @param R2 correlation/covariance/rank-based correlation matrix of dataset \code{X2}.
#' @param R12 correlation/covariance/rank-based correlation matrix between dataset \code{X1} and dataset \code{X2}.
#' @param lambda1 tuning parameter dataset \code{X1}. a scalar value.
#' @param lambda2 tuning parameter dataset \code{X2}. a scalar value.
#' @param tol tolerance for eigenvalues. Refer to \code{standardCCA} function.
#'
#' @return \code{myrcc} returns a data.frame containing
#' \itemize{
#'       \item{cancor: }{estimated canonical correlation.}
#'       \item{w1: }{estimated canonical direction \eqn{w1}.}
#'       \item{w2: }{estimated canonical direction \eqn{w2}.}
#' }
#' @export
#'
myrcc <- function(R1, R2, R12, lambda1, lambda2, tol = 1e-4){

  C1 <- R1 + diag(lambda1, ncol(R1))
  C2 <- R2 + diag(lambda2, ncol(R2))
  C12 <- R12
  res <- standardCCA(S1 = C1, S2 = C2, S12 = C12, tol = tol)

  return(list(cancor = res$cancor,
              w1 = res$w1,
              w2 = res$w2))
}

