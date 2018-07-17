#' Construct a correlation matrix
#'
#' Functions to create autocorrelation matrix of size p with parameter rho and block correlation matrix (p by p) using group index (of length p) and different parameter rho for each group.
#' @name corrmat
#' @rdname corrmat
#' @aliases autocor
#' @param p Dimension of matrix
#' @param rho Correlation
#' @export
autocor <- function(p, rho){
  Sigma <- rho^abs(outer(1:p, 1:p, "-"))
  return(Sigma)
}
#' @rdname corrmat
#' @aliases blockcor
#' @inheritParams rho
#' @param groupind Group index indicating belonging to which block.
#' @seealso See \code{\link{GenerateData}} for an example.
#' @export
blockcor <- function(rho, groupind){
  # groupind variable is a vector of indices. For example, c(rep(1,3), c(2,5)) for two groups. p=8.
  # if (p!=length(group)){
  #   stop("p and groupind must have the same length.")
  # }
  p <- length(groupind)
  blk <- unique(groupind)
  if (length(rho) != length(blk)){
    if(length(rho) == 1){
      rho <- rep(rho, length(blk))
    } else {
      stop("rho and kinds of group must match.")
    }
  }
  Sigma <- matrix(0, p, p)

  for (j in 1:length(blk)){
    coef <- which(groupind %in% blk[j])
    Sigma[coef, coef] <- rho[j]
  }
  diag(Sigma) = 1
  return(Sigma)
}

#' Mixed type simulation data generator for sparse CCA
#'
#' \code{GenerateData} is used to generate two sets of data of mixed types for sparse CCA under the Gaussian copula model.
#'
#' @param n Sample size
#' @param trueidx1 True canonical direction of length p1 for \code{X1}. It will be automatically normalized such that \eqn{w_1^T \Sigma_1 w_1 = 1}.
#' @param trueidx2 True canonical direction of length p2 for \code{X2}. It will be automatically normalized such that \eqn{w_2^T \Sigma_2 w_2 = 1}.
#' @param Sigma1 True correlation matrix of latent variable \code{Z1} (p1 by p1).
#' @param Sigma2 True correlation matrix of latent variable \code{Z2} (p2 by p2).
#' @param maxcancor True canonical correlation between \code{Z1} and \code{Z2}.
#' @param copula1 Copula type for first dataset. U1 = f(Z1). Currently, we have two options: "exp", "cube".
#' @param copula2 Copula type for second dataset. U2 = f(Z2).
#' @param type1 Type of first dataset \code{X1}.
#' @param type2 Type of second dataset \code{X2}.
#' @param muZ Mean of multivariate normal for latent data generation.
#' @param c1 Constant threshold for \code{X1} only needed for "trunc" and "binary" data type.
#' @param c2 Constant threshold for \code{X2} only needed for "trunc" and "binary" data type.
#'
#' @return \code{GenerateData} returns a data.frame containing
#' \itemize{
#'       \item{Z1: }{latent numeric data matrix (n by p1).}
#'       \item{Z2: }{latent numeric data matrix (n by p2).}
#'       \item{X1: }{observed numeric data matrix (n by p1).}
#'       \item{X2: }{observed numeric data matrix (n by p2).}
#'       \item{true_w1: }{normalized true canonical direction of length p1 for \code{X1}.}
#'       \item{true_w2: }{normalized true canonical direction of length p2 for \code{X2}.}
#'       \item{type: }{a vector containing types of two datasets.}
#'       \item{maxcancor: }{true canonical correlation between \code{Z1} and \code{Z2}.}
#'       \item{c1: }{constant threshold for \code{X1} for "trunc" and "binary" data type.}
#'       \item{c2: }{constant threshold for \code{X2} for "trunc" and "binary" data type.}
#'       \item{Sigma: }{true latent correlation matrix of \code{Z1} and \code{Z2} ((p1+p2) by (p1+p2)).}
#' }
#' @export
#'
#' @importFrom MASS mvrnorm
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
GenerateData <- function(n, trueidx1, trueidx2, Sigma1, Sigma2, maxcancor,
                         copula1 = "no", copula2 = "no",
                         type1 = "continuous", type2 = "continuous", muZ = NULL, c1 = NULL, c2 = NULL
                         # , seed = 1
){

  type <- rep(NA, 2); type[1] <- type1; type[2] <- type2
  p1 <- length(trueidx1)
  p2 <- length(trueidx2)
  p <- p1 + p2

  # normalize to satisfy t(theta1)%*%Sigma1%*%theta1=1
  th1 <- trueidx1/sqrt(as.numeric(crossprod(trueidx1, Sigma1 %*% trueidx1)))
  th2 <- trueidx2/sqrt(as.numeric(crossprod(trueidx2, Sigma2 %*% trueidx2)))
  Sigma12 <- maxcancor*Sigma1%*%th1%*%t(th2)%*%Sigma2
  JSigma <- rbind(cbind(Sigma1, Sigma12), cbind(t(Sigma12), Sigma2))

  # jointly generate X and Y using two canonical pairs
  # set.seed(seed)
  if (is.null(muZ)) { muZ <- rep(0, p) }
  dat <- MASS::mvrnorm(n, mu = muZ, Sigma = JSigma)

  Z1 <- dat[, 1:p1]
  Z2 <- dat[, (p1+1):p]

  if(copula1 != "no"){
    if(copula1 == "exp"){
      Z1 <- exp(Z1)
    }else if(copula1 == "cube"){
      Z1 <- Z1^3
    }
  }
  if(copula2 != "no"){
    if(copula2 == "exp"){
      Z2 <- exp(Z2)
    }else if(copula2 == "cube"){
      Z2 <- Z2^3
    }
  }

  if(type1 != "continuous"){
    if(length(c1) != p1) { warning("The length of threshold vector does not match with the size of the data.") }
    if(length(c1) == 1) { warning("Same therhold is applied to the all variables in the first set.") }
  }
  if(type2 != "continuous"){
    if(length(c2) != p2) { warning("The length of threshold vector does not match with the size of the data.") }
    if(length(c2) == 1) { warning("Same therhold is applied to the all variables in the second set.") }
  }

  if(type1 == "continuous") { X1 <- Z1 }else if(type1 == "trunc") { X1 <- ifelse(Z1 > c1, Z1, 0) }else if (type1 == "binary") { X1 <- ifelse(Z1 > c1, 1, 0) }
  if(type2 == "continuous") { X2 <- Z2 }else if(type2 == "trunc") { X2 <- ifelse(Z2 > c2, Z2, 0) }else if (type2 == "binary") { X2 <- ifelse(Z2 > c2, 1, 0) }

  return(list(Z1 = Z1, Z2 = Z2, X1 = X1, X2 = X2, true_w1 = th1, true_w2 = th2, type = type, maxcancor = maxcancor, c1 = c1, c2 = c2, Sigma = JSigma))
}


