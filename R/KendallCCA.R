lambdaseq_generate <- function(nlambda = 20, initlam1 = NULL, initlam2 = NULL, lam.eps = 1e-02, Sigma1, Sigma2, Sigma12, w1init = NULL, w2init = NULL){
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


#' Sparse CCA for data of mixed types with BIC criterion
#'
#' Applies sparse canonical correlation analysis (CCA) for high-dimensional data of mixed types (continuous/biary/truncated continuous). Derived rank-based estimator instead of sample correlation matrix is implemented. There are two types of BIC criteria for variable selection. We found that BIC1 works best for variable selection, whereas BIC2 works best for prediction.
#'
#' @param X1 A numeric data matrix (n by p1).
#' @param X2 A numeric data matrix (n by p2).
#' @param type1 A type of data \code{X1} among "continuous", "binary", "trunc".
#' @param type2 A type of data \code{X2} among "continuous", "binary", "trunc".
#' @param lamseq1 A tuning parameter sequence for \code{X1}. The length should be the same as \code{lamseq2}.
#' @param lamseq2 A tuning parameter sequence for \code{X2}. The length should be the same as \code{lamseq1}.
#' @param initlam1 An initial value to generate a tuning parameter sequence for \code{X1}, which is the maximum value of tuning parameter grid.
#' @param initlam2 An initial value to generate a tuning parameter sequence for \code{X2}, which is the maximum value of tuning parameter grid.
#' @param nlambda The number of tuning parameter sequence lambda - default is 20.
#' @param lam.eps A ratio of the smallest value for lambda to the maximum value of lambda.
#' @param w1init An initial vector of length p1 for canonical direction \eqn{w1}.
#' @param w2init An initial vector of length p2 for canonical direction \eqn{w2}.
#' @param BICtype Either 1 or 2: For more details for two options, see the reference.
#' @param KendallR An estimated Kendall \eqn{\tau} matrix. The default is NULL, which means that it will be automatically estimated by Kendall's \eqn{\tau} estimator unless the user supplies.
#' @param tol The desired accuracy (convergence tolerance).
#' @param maxiter The maximum number of iterations allowed.
#' @param verbose If \code{verbose = FALSE}, printing convergence error when the convergence is failed after \code{maxiter} is disabled. The default value is \code{TRUE}.
#'
#' @references
#' Yoon G., Carroll R.J. and Gaynanova I. (2018+) \href{http://arxiv.org/abs/1807.05274}{"Sparse semiparametric canonical correlation analysis for data of mixed types"}.
#' @return \code{mixedCCA} returns a data.frame containing
#' \itemize{
#'       \item{KendallR: }{estimated Kendall's \eqn{\tau} matrix estimator.}
#'       \item{lambda_seq: }{the values of \code{lambda} used for sparse CCA.}
#'       \item{w1: }{estimated canonical direction \eqn{w1}.}
#'       \item{w2: }{estimated canonical direction \eqn{w2}.}
#'       \item{cancor: }{estimated canonical correlation.}
#'       \item{selected_x1: }{indices of selected variables in \code{X1}.}
#'       \item{selected_x2: }{indices of selected variables in \code{X2}.}
#' }
#' @importFrom Rcpp sourceCpp
#' @import RcppArmadillo
#' @import stats
#'
#' @example man/examples/mixedCCA_ex.R
#' @useDynLib mixedCCA
#' @export
mixedCCA <- function(X1, X2, type1, type2, lamseq1 = NULL, lamseq2 = NULL, initlam1 = NULL, initlam2 = NULL,
                     nlambda = 20, lam.eps = 1e-02,
                     w1init = NULL, w2init = NULL, BICtype,
                     KendallR = NULL,
                     tol = 1e-3, maxiter = 1000, verbose = TRUE){
  n <- nrow(X1)
  p1 <- ncol(X1); p2 <- ncol(X2);
  p <- p1 + p2

  ### Compute Kendall tau.
  if(is.null(KendallR)){
    R <- estimateR_mixed(X1, X2, type1 = type1, type2 = type2)
  } else { R <- KendallR; rm(KendallR) }

  R1 <- as.matrix(R$R1)
  R2 <- as.matrix(R$R2)
  R12 <- as.matrix(R$R12)

  if (is.null(w1init) | is.null(w2init)){
      res.regul <- estim.regul_crossvalidation(X1, X2)
      RCCA <- myrcc(X1, X2, res.regul$lambda1.optim, res.regul$lambda2.optim)

      if (is.null(w1init)){ w1init <- as.matrix(RCCA$w1, ncol=1) }
      if (is.null(w2init)){ w2init <- as.matrix(RCCA$w2, ncol=1) }
  }

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
                           w1init, w2init, BICtype = BICtype, verbose = verbose)
    if(sum(coefall[1:p1]!=0)==0){
      lamtmp[[1]] <- exp(seq(log(lamtmp[[1]][2]), log(lambda_seq[[1]][nlambda]), length.out = nlambda))
    }
    if(sum(coefall[(p1+1):(p1+p2)]!=0)==0){
      lamtmp[[2]] <- exp(seq(log(lamtmp[[2]][2]), log(lambda_seq[[2]][nlambda]), length.out = nlambda))
    }
  }
  lambda_seq <- lamtmp

  coeff <- find_w12bic(n, R1, R2, R12, lambda_seq[[1]], lambda_seq[[2]], maxiter = maxiter, tol = tol,
                       w1init, w2init, BICtype = BICtype, verbose = verbose)

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


# CV for canonical ridge method
estim.regul_crossvalidation <- function (X1, X2, lambda1grid = NULL, lambda2grid = NULL, nfolds = 5){

  if (is.null(lambda1grid)) {
    lambda1grid <- matrix(seq(0.001, 1, length = 5),nrow=1)
  } else {
    lambda1grid <- matrix(lambda1grid, nrow=1)
  }
  if (is.null(lambda2grid)) {
    lambda2grid <- matrix(seq(0.001, 1, length = 5),nrow=1)
  } else {
    lambda2grid <- matrix(lambda2grid, nrow=1)
  }

  lambda1.matrix <- matrix(rep(lambda1grid,length(lambda1grid)), ncol=length(lambda2grid), byrow=T)
  lambda2.matrix <- matrix(sort(rep(lambda2grid,length(lambda2grid))), ncol=length(lambda1grid), byrow=T)

  cvscores <- apply(lambda1grid, 2, l1function, X1=X1, X2=X2, lambda2grid=lambda2grid, nfolds=nfolds) #cv-score
  cv.optim <- cvscores[which.max(cvscores)]
  lambda1.optim <- lambda1.matrix[which.max(cvscores)]
  lambda2.optim <- lambda2.matrix[which.max(cvscores)]

  ##OUTPUT
  out = list(lambda1.optim=lambda1.optim, lambda2.optim=lambda1.optim, cv.optim=cv.optim)
  return(out)
}

l1function<-function(lam1, X1, X2, lambda2grid, nfolds){ # AUXILIARY FUNCTION CANONICAL RIDGE
  testcor <- apply(lambda2grid, 2, l2function, X1=X1, X2=X2, lambda1fixed=lam1, nfolds=nfolds)
  return(testcor)
}

l2function <- function(lam2, X1, X2, lambda1fixed, nfolds){ # AUXILIARY FUNCTION CANONICAL RIDGE
  RCCcvm <- RCC_crossvalidation(X1=X1, X2=X2, lambda1=lambda1fixed, lambda2=lam2, nfolds=nfolds)
  return(RCCcvm)
}


RCC_crossvalidation <- function (X1, X2, lambda1, lambda2, nfolds) { # AUXILIARY FUNCTION CANONICAL RIDGE: n.cv-fold cross-validation

  n <- nrow(X1)
  id <- sample(rep(seq_len(nfolds), length.out = n))

  cv <- 0
  # For each fold, do
  for (nf in 1: nfolds){
    ### set training data and test data
    xtrain1 <- X1[id != nf, ]
    xtrain2 <- X2[id != nf, ]
    xtest1 <- X1[id == nf, ]
    xtest2 <- X2[id == nf, ]

    ### compute Canonical Ridge
    res = myrcc(xtrain1, xtrain2, lambda1, lambda2)

    #### Calculate cv metrics for each lambda
    xscore = xtest1 %*% res$w1
    yscore = xtest2 %*% res$w2
    cv <- cv + abs(cor(xscore,yscore,use="pairwise"))
  }
  cvm <- sum(cv)/nfolds
  return(cvm)

}
