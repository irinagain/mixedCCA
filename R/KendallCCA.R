# sparse CCA code using Rcpp.


#' Generate tuning parameter sequence for BIC.
#'
#' @param nlambda the number of tuning parameter sequence \code{lamseq} - default is 100.
#' @param initlam1 initial value to generate a tuning parameter sequence for \eqn{X1}.
#' @param initlam2 initial value to generate a tuning parameter sequence for \eqn{X2}.
#' @param lam.eps ratio of the smallest value for lambda to the maximum value of lambda.
#' @param Sigma1 correlation matrix of \eqn{X1}.
#' @param Sigma2 correlation matrix of \eqn{X2}.
#' @param Sigma12 correlation matrix between \eqn{X1} and \eqn{X2}.
#' @param w1init initial value of canonical direction \eqn{w1}, which helps estimate \code{initlam1} better. The default is same coordinates for all elements such that \eqn{w^T \Sigma w = 1}.
#' @param w2init initial value of canonical direction \eqn{w2}, which helps estimate \code{initlam2} better. The default is same coordinates for all elements such that \eqn{w^T \Sigma w = 1}.
#'
#' @return a list containing two vectors of lambda sequences for \eqn{X1} and \eqn{X2}.
#'
lambdaseq_generate <- function(nlambda = 100, initlam1 = NULL, initlam2 = NULL, lam.eps = 1e-02, Sigma1, Sigma2, Sigma12, w1init = NULL, w2init = NULL){
  init.lambda <- rep(NA, 2)
  p1 <- nrow(Sigma1); p2 <- nrow(Sigma2)
  if (is.null(initlam1)) {
    if (is.null(w2init)) {
      w2init <- rep(1,p2)/sqrt(as.numeric(crossprod(rep(1,p2), Sigma2 %*% rep(1,p2))))
    }
    init.lambda[1] <- max(abs(Sigma12%*%w2init)) # estimated w1 for different (but close) lambda could be used
  }

  if (is.null(initlam2)) {
    if (is.null(w1init)) {
      w1init <- rep(1,p1)/sqrt(as.numeric(crossprod(rep(1,p1), Sigma1 %*% rep(1,p1))))
    }
    init.lambda[2] <- max(abs(t(Sigma12)%*%w1init)) # estimated w2 for different (but close) lambda could be used
  }

  lambda <- list("vector")
  # lambda values in decreasing order.
  lambda[[1]] <- exp(seq(log(init.lambda[1]), log(lam.eps*init.lambda[1]), length.out = nlambda))
  lambda[[2]] <- exp(seq(log(init.lambda[2]), log(lam.eps*init.lambda[2]), length.out = nlambda))

  return(lambda)
}


#' Kendall sparse CCA with BIC criterion
#'
#' @param X1 numeric data matrix (n by p1).
#' @param X2 numeric data matrix (n by p2).
#' @param type1 the type of data matrix \code{X1} among "continuous", "binary", "trunc".
#' @param type2 the type of data matrix \code{X2} among "continuous", "binary", "trunc".
#' @param lamseq1 tuning parameter sequence for \code{X1}. The length should be the same as lamseq2.
#' @param lamseq2 tuning parameter sequence for \code{X2}. The length should be the same as lamseq1.
#' @param initlam1 initial value to generate a tuning parameter sequence for \code{X1}.
#' @param initlam2 initial value to generate a tuning parameter sequence for \code{X2}.
#' @param nlambda the number of tuning parameter sequence \code{lamseq} - default is 100.
#' @param lam.eps ratio of the smallest value for lambda to the maximum value of lambda.
#' @param w1init initial value of canonical direction \eqn{w1}.
#' @param w2init initial value of canonical direction \eqn{w2}.
#' @param BICtype 1 or 2: For more details for two options, see the reference.
#' @param KendallR estimated Kendall \eqn{\tau} matrix. The default is NULL, which means that it will be automatically estimated by Kendall's \eqn{\tau} estimator unless the user supplies.
#' @param tol the desired accuracy (convergence tolerance).
#' @param maxiter maximum number of iterations allowed.
#'
#' @references Yoon G, Carroll R and Gaynanova I. Sparse semiparametric canonical correlation analysis for high-dimensional data of mixed types. (2018+)
#' @return A data.frame containing
#' \itemize{
#'       \item{KendallR: }{estimated Kendall's \eqn{\tau} matrix estimator.}
#'       \item{lambda_seq: }{the values of \code{lambda} used for sparse CCA.}
#'       \item{w1: }{estimated canonical direction \eqn{w1}.}
#'       \item{w2: }{estimated canonical direction \eqn{w2}.}
#'       \item{cancor: }{estimated canonical correlation.}
#'       \item{selected_x1: }{indices of selected variables in \code{X1}.}
#'       \item{selected_x2: }{indices of selected variables in \code{X2}.}
#' }
#' @useDynLib mixedCCA
#' @importFrom Rcpp sourceCpp
#' @import RcppArmadillo
#'
#' @examples
#' n <- 100; p1 <- 5; p2 <- 7
#' maxcancor <- 0.9
#' Sigma1 <- autocor(p1, 0.7)
#' groupind <- c(rep(1, 2), rep(2, p2-2))
#' Sigma2 <- blockcor(0.7, groupind)
#' copula1 <- "exp"; copula2 <- "cube"
#' mu <- rep(0, p1+p2)
#' type1 <- type2 <- "trunc"
#' c1 <- rep(0, p1); c2 <- rep(0, p2)
#' trueidx1 <- c(0, 0, 1, 1, 1)
#' trueidx2 <- c(1, 0, 1, 0, 0, 1, 1)
#' simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2, maxcancor = maxcancor,
#'                        Sigma1 = Sigma1, Sigma2 = Sigma2,
#'                        copula1 = copula1, copula2 = copula2,
#'                        muZ = mu,
#'                        type1 = type1, type2 = type2, c1 = c1, c2 = c2
#' )
#'
#' X1 <- simdata$X1
#' X2 <- simdata$X2
#'
#' RCCA <- myrcc(X1, X2, 0.01, 0.1)
#' w1init <- RCCA$w1[, 1]
#' w2init <- RCCA$w2[, 1]

#' # Kendall CCA with BIC1
#' kendallcca1 <- mixedCCA(X1, X2, type1 = type1, type2 = type2, lam.eps = 0.01, nlambda = 20,
#'                        w1init = w1init, w2init = w2init, BICtype = 1)
#'
#' @export
mixedCCA <- function(X1, X2, type1, type2, lamseq1 = NULL, lamseq2 = NULL, initlam1 = NULL, initlam2 = NULL, nlambda = 100, lam.eps = 1e-02,
                    w1init, w2init, BICtype,
                    KendallR = NULL,
                    tol = 1e-3, maxiter = 1000){
  n <- nrow(X1)
  p1 <- ncol(X1); p2 <- ncol(X2);
  p <- p1 + p2

  ### Compute Kendall tau.
  if(is.null(KendallR)){
    R <- estimateR_mixed(X1, X2, type1 = type1, type2 = type2)
    cat("Complete the Kendall tau calculation.\n")
  } else { R <- KendallR; rm(KendallR) }

  R1 <- as.matrix(R$R1)
  R2 <- as.matrix(R$R2)
  R12 <- as.matrix(R$R12)

  ### Create lambda sequences
  lambda_seq <- list()
  if (is.null(lamseq1) & is.null(lamseq2)){
    lambda_seq <- lambdaseq_generate(nlambda = nlambda, initlam1 = initlam1, initlam2 = initlam2, lam.eps = lam.eps,
                                     Sigma1 = R$R1, Sigma2 = R$R2, Sigma12 = R$R12, w1init = w1init, w2init = w2init)
  } else {
    if(length(lamseq1) == length(lamseq2)){
    nlambda <- length(lamseq1)
    lambda_seq[[1]] <- lamseq1
    lambda_seq[[2]] <- lamseq2
    } else { warning( "The lengths of lambda sequences for two variables are different." ); stop; }
  }

  ### To remove lambda values which make all elements have zero coefficients.
  coefall <- rep(0, p)
  lamtmp <- lambda_seq
  while(sum(coefall[1:p1]!=0)==0 | sum(coefall[(p1+1):(p1+p2)]!=0)==0){
    coefall <- find_w12bic(n, R1, R2, R12, rep(lamtmp[[1]][1], p1), rep(lamtmp[[2]][1], p2), maxiter = maxiter, tol = tol,
                           w1init, w2init, BICtype = BICtype)
    if(sum(coefall[1:p1]!=0)==0){
      lamtmp[[1]] <- exp(seq(log(lamtmp[[1]][2]), log(lambda_seq[[1]][nlambda]), length.out = nlambda))
    }
    if(sum(coefall[(p1+1):(p1+p2)]!=0)==0){
      lamtmp[[2]] <- exp(seq(log(lamtmp[[2]][2]), log(lambda_seq[[2]][nlambda]), length.out = nlambda))
    }
  }
  lambda_seq <- lamtmp

  coeff <- find_w12bic(n, R1, R2, R12, lambda_seq[[1]], lambda_seq[[2]], maxiter = maxiter, tol = tol,
                       w1init, w2init, BICtype = BICtype)

  w1 <- coeff[1:p1]
  w2 <- coeff[(p1+1):(p1+p2)]
  ix1 <- which(w1!=0)
  ix2 <- which(w2!=0)

  cancor <- t(w1)%*%R12%*%w2

  return(list(KendallR = R,
              lambda_seq = lambda_seq,
              w1 = w1, w2 = w2,
              cancor = cancor,
              selected_x1 = ix1, selected_x2 = ix2))
}


#' Classical CCA using correlation matrix estimator.
#'
#' @param S1 correlation matrix of \eqn{X1}.
#' @param S2 correlation matrix of \eqn{X2}.
#' @param S12 correlation matrix between \eqn{X1} and \eqn{X2}.
#' @param tol cutoff for eigenvalues. During the calculation, eigenvalues smaller than cutoff are set to zero.
#'
#' @return A data.frame containing
#' \itemize{
#'       \item{cancor: }{estimated canonical correlation.}
#'       \item{w1: }{estimated canonical direction \eqn{w1}.}
#'       \item{w2: }{estimated canonical direction \eqn{w2}.}
#' }
#' @export
#'
#' @importFrom irlba irlba
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

  outsvd <- irlba(S1isqrt %*% S12 %*% S2isqrt, nv = 1) # svd2(S1isqrt %*% S12 %*% S2isqrt, nu = 1, nv = 1)

  # We are considering only one pair of canonical covariates.
  cancor <- max(outsvd$d)

  w1 <- S1isqrt %*% outsvd$u
  w2 <- S2isqrt %*% outsvd$v

  return(list(cancor = cancor, w1 = w1, w2 = w2))
}


#' Canonical Ridge
#' This code is modified from rcc in the package \code{CCA}.
#'
#' @param X1 numeric data matrix (n by p1).
#' @param X2 numeric data matrix (n by p2).
#' @param lambda1 regularization parameter sequence for \code{X1}. The length should be the same as \code{ncol(X1)}.
#' @param lambda2 regularization parameter sequence for \code{X2}. The length should be the same as \code{ncol(X2)}.
#' @return A data.frame containing
#' \itemize{
#'       \item{cancor: }{estimated canonical correlation.}
#'       \item{w1: }{estimated canonical direction \eqn{w1}.}
#'       \item{w2: }{estimated canonical direction \eqn{w2}.}
#' }
#' @export
#'
#' @seealso \code{\link{mixedCCA}}
myrcc <- function(X1, X2, lambda1, lambda2){
  X1names <- dimnames(X1)[[2]]
  X2names <- dimnames(X2)[[2]]
  ind.names <- dimnames(X1)[[1]]
  C1 <- var(X1, na.rm = TRUE, use = "pairwise") + diag(lambda1, ncol(X1))
  C2 <- var(X2, na.rm = TRUE, use = "pairwise") + diag(lambda2, ncol(X2))
  C12 <- cov(X1, X2, use = "pairwise")
  res <- standardCCA(C1, C2, C12)

  return(list(cancor = res$cancor,
              w1 = res$w1,
              w2 = res$w2))
}


