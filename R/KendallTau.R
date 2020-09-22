#' Kendall's tau correlation
#'
#' Calculate Kendall's tau correlation.
#' \deqn{ \hat{\tau}_{jk} = \frac{2}{n(n-1)}\sum_{1\le i<i'\le n} sign(X_{ji}-X_{ji'}) sign(X_{ki}-X_{ki'}) }
#' The function \code{KendallTau} calculates Kendall's tau correlation between two variables, returning a single correlation value. The function \code{Kendall_matrix} returns a correlation matrix.
#'
#' @param x A numeric vector.
#' @param y A numeric vector.
#'
#' @rdname KendallTau
#' @return \code{KendallTau(x, y)} returns one Kendall's tau correlation value between two vectors, \code{x} and \code{y}.
#' @examples
#'
#' n <- 100 # sample size
#' r <- 0.8 # true correlation
#'
#' ### vector input
#' # Data generation (X1: truncated continuous, X2: continuous)
#' Z <- mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, r, r, 1), nrow = 2))
#' X1 <- Z[,1]
#' X1[Z[,1] < 1] <- 0
#' X2 <- Z[,2]
#'
#' KendallTau(X1, X2)
#' Kendall_matrix(X1, X2)
#'
#' ### matrix data input
#' p1 <- 3; p2 <- 4 # dimension of X1 and X2
#' JSigma <- matrix(r, nrow = p1+p2, ncol = p1+p2); diag(JSigma) <- 1
#' Z <- mvrnorm(n, mu = rep(0, p1+p2), Sigma = JSigma)
#' X1 <- Z[,1:p1]
#' X1[Z[,1:p1] < 0] <- 0
#' X2 <- Z[,(p1+1):(p1+p2)]
#'
#' Kendall_matrix(X1, X2)
#'
#' @importFrom pcaPP cor.fk
#' @export
KendallTau <- function(x, y){ # both x and y are vectors, not matrix.
  # Based on cor.fk function from pcaPP package to make the computation faster.
  # It can handle ties.
  if (length(x) != length(y)){ # Check of they have the same length
    stop ("x and y must have same length.")
  }
  if (sum(is.na(x)) + sum(is.na(y)) > 0){
    index_new = which(!is.na(x) & !is.na(y))
    if (length(index_new) < 2){
      stop("Too many missing measuments in one of the variables, filter for missing data.")
    }else{
      x <- x[index_new]
      y <- y[index_new]
    }
  }

  n <- length(x)
  n0 <- n * (n - 1) / 2
  if (length(unique(x)) != n) {
    x <- as.vector(x) # sometimes input x is a matrix n by 1, which gives errors for rle function below.
    x.info <- rle(sort(x))
    t1 <- x.info$lengths[x.info$lengths > 1]
    n1 <- sum(t1 * (t1 - 1) / 2)
  } else {
    n1 <- 0
  }
  if (length(unique(y)) != n) {
    y <- as.vector(y) # sometimes input y is a matrix n by 1, which gives errors for rle function below.
    y.info <- rle(sort(y))
    u1 <- y.info$lengths[y.info$lengths > 1]
    n2 <- sum(u1 * (u1 - 1) / 2)
  } else {
    n2 <- 0
  }

  tau <- pcaPP::cor.fk(x, y) * sqrt(n0 - n1) * sqrt(n0 - n2) / n0

  return(tau)
}

#' @param X A numeric matrix (n by p1).
#' @param Y A numeric matrix (n by p2).
#'
#' @rdname KendallTau
#'
#' @return \code{Kendall_matrix(X)} returns a p1 by p1 matrix of Kendall's tau correlation coefficients. \code{Kendall_matrix(X, Y)} returns a p1 by p2 matrix of Kendall's tau correlation coefficients.
#'
#' @export
#' @importFrom pcaPP cor.fk
Kendall_matrix <- function(X, Y = NULL){ # X and Y are matrix.
  if (is.null(Y)) {
    X <- as.matrix(X) # In case that X is vector, this line gives you one column matrix.
    p1 <- ncol(X)
    tau <- matrix(1, p1, p1)
    if (p1 > 1){
      for (i in 1:(p1 - 1)){
        for (j in (i + 1):p1){ # calculate KendallTau elementwise.
          tau[i, j] <- tau[j, i] <- KendallTau(X[, i], X[, j])
        }
      }
    }
  } else {
    X <- as.matrix(X); Y <- as.matrix(Y) # In case that X and Y are vectors, this line gives you one column matrix.
    p1 <- ncol(X); p2 <- ncol(Y)
    if (p1 == 1 & p2 == 1){
      tau <- KendallTau(X, Y)
    } else {
      tau <- matrix(0, p1, p2)
      for (i in 1:p1){
        for (j in 1:p2){ # calculate KendallTau elementwise.
          tau[i, j] <- KendallTau(X[, i], Y[, j])
        }
      }
    }
  }
  return(tau)
}
