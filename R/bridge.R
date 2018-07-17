# bridge function for mixed data depending on the type

#' Bridge functions
#' @aliases bridge
#'
#'
#' Bridge function \eqn{F} connects Kendall's tau \eqn{\tau} estimation to unbiased latent correlation matrix \eqn{\Sigma} under Gaussian copula model. There are three possible types of variables: continuos, binary, truncated continuous. Depending on the combination of different types of variables, different bridge function should be applied. \code{bridge_select} shows which bridge function should be used based on data types. Taking the first letter from the type of the variables, the combination of the variables' type follows the bridgeF. For example, bc indicates the first dataset is binary and the second dataset is continuous.
#'
#' \deqn{ F( \Sigma_{jk} ) = \tau_{jk} }
#' The input \eqn{r} is an element of \eqn{\Sigma} and the bridge function is calculated elementwise.
#'
#' @param r correlation.
#' @param type1 Type of the dataset \eqn{X_{1}}. Choose among "continous", "binary", "trunc".
#' @param type2 Type of the dataset \eqn{X_{2}}.
#' @param zratio1 proportion of zero values of dataset \eqn{X_{1}}.
#' @param zratio2 proportion of zero values of dataset \eqn{X_{2}}.
#'
#' @rdname bridge
#' @examples
#'
#' set.seed(1)
#' n <- 100 # sample size
#' r <- 0.8 # true correlation
#'
#' # Data generation
#'
#' # X1: truncated continuous, X2: continuous
#' Z <- mvrnorm(n, mu = c(0, 0), Sigma = matrix(c(1, r, r, 1), nrow = 2))
#' X1 <- Z[,1]
#' X1[Z[,1] < 1] <- 0
#' X2 <- Z[,2]
#' zratio1 <- mean(X1==0)
#' # Choose bridge function
#' bridge <- bridge_select(type1 = "trunc", type2 = "continuous")
#' K12 <- KendallTau(X1, X2)
#' K12
#' f1 <- function(r)(bridge(r, zratio1 = zratio1) - K12)^2
#' optimize(f1, lower = -0.99, upper = 0.99, tol = 0.001)[[1]]
#'
#' # X1: truncated continuous, X2: truncated continuous
#' X2[Z[,2] < 1] <- 0
#' zratio2 <- mean(X2==0)
#' bridge <- bridge_select(type1 = "trunc", type2 = "trunc")
#' K12 <- KendallTau(X1, X2)
#' K12
#' f1 <- function(r)(bridge(r, zratio1 = zratio1, zratio2 = zratio2) - K12)^2
#' optimize(f1, lower = -0.99, upper = 0.99, tol = 0.001)[[1]]
#'
#' # X1: truncated continuous, X2: binary
#' X2[X2 != 0] <- 1
#' zratio2 <- mean(X2==0)
#' bridge <- bridge_select(type1 = "trunc", type2 = "binary")
#' K12 <- KendallTau(X1, X2)
#' K12
#' f1 <- function(r)(bridge(r, zratio1 = zratio1, zratio2 = zratio2) - K12)^2
#' optimize(f1, lower = -0.99, upper = 0.99, tol = 0.001)[[1]]
#' @export
bridge_select <- function(type1 = "trunc", type2 = "continuous") {
  if (type1 == "binary" & type2 == "binary") { bridge_select <- bridgeF_bb
  } else if (type1 == "trunc" & type2 == "trunc") { bridge_select <- bridgeF_tt
  } else if (type1 == "trunc" & type2 == "continuous") { bridge_select <- bridgeF_tc
  } else if (type1 == "continuous" & type2 == "trunc") { bridge_select <- bridgeF_ct
  } else if (type1 == "binary" & type2 == "continuous") { bridge_select <- bridgeF_bc
  } else if (type1 == "continuous" & type2 == "binary") { bridge_select <- bridgeF_cb
  } else if (type1 == "trunc" & type2 == "binary") { bridge_select <- bridgeF_tb
  } else if (type1 == "binary" & type2 == "trunc") { bridge_select <- bridgeF_bt
  } else {
    stop("Unrecognized type of variables. Should be one of continuous, binary or trunc.")
  }
}
#' @rdname bridge
#' @import stats
#' @importFrom fMultivar pnorm2d
#' @importFrom mnormt pmnorm
bridgeF_bc <- function(r, zratio1, zratio2 = NULL){
  # binary and continuous
  de1 <- qnorm(zratio1)
  res <- as.numeric( 4*fMultivar::pnorm2d(de1, 0, rho = r/sqrt(2)) - 2*zratio1 )
  return(res)
}
#' @rdname bridge
bridgeF_cb <- function(r, zratio1 = NULL, zratio2){
  # continuous and binary
  de2 <- qnorm(zratio2)
  res <- as.numeric( 4*fMultivar::pnorm2d(0, de2, rho = r/sqrt(2)) - 2*zratio2 )
  return(res)
}
#' @rdname bridge
#' @import stats
#' @importFrom fMultivar pnorm2d
#' @importFrom mnormt pmnorm
bridgeF_tb <- function(r, zratio1, zratio2){
  # truncated and binary
  de1 <- qnorm(zratio1)
  de2 <- qnorm(zratio2)
  mat1 <- matrix(c(1, -r, 1/sqrt(2),
                   -r, 1, -r/sqrt(2),
                   1/sqrt(2), -r/sqrt(2), 1), nrow = 3)
  mat2 <- matrix(c(1, 0, -1/sqrt(2),
                   0, 1, -r/sqrt(2),
                   -1/sqrt(2), -r/sqrt(2), 1), nrow = 3)
  res <- as.numeric(
    2*(1-zratio1)*(zratio2)-
      2*mnormt::pmnorm(c(-de1, de2, 0), mean = rep(0, 3), varcov = mat1)-
      2*mnormt::pmnorm(c(-de1, de2, 0), mean = rep(0, 3), varcov = mat2)
  )
  return(res)
}
#' @rdname bridge
#' @import stats
#' @importFrom fMultivar pnorm2d
#' @importFrom mnormt pmnorm
bridgeF_bt <- function(r, zratio1, zratio2){
  # binary and truncated
  de1 <- qnorm(zratio2)
  de2 <- qnorm(zratio1)
  mat1 <- matrix(c(1, -r, 1/sqrt(2),
                   -r, 1, -r/sqrt(2),
                   1/sqrt(2), -r/sqrt(2), 1), nrow = 3)
  mat2 <- matrix(c(1, 0, -1/sqrt(2),
                   0, 1, -r/sqrt(2),
                   -1/sqrt(2), -r/sqrt(2), 1), nrow = 3)
  res <- as.numeric(
    2*(1-zratio2)*(zratio1)-
      2*mnormt::pmnorm(c(-de1, de2, 0), mean = rep(0, 3), varcov = mat1)-
      2*mnormt::pmnorm(c(-de1, de2, 0), mean = rep(0, 3), varcov = mat2)
  )
  return(res)
}
#' @rdname bridge
#' @import stats
#' @importFrom fMultivar pnorm2d
#' @importFrom mnormt pmnorm
bridgeF_tc <- function(r, zratio1, zratio2 = NULL){
  # truncated and continuous
  de1 <- qnorm(zratio1)
  mat2 <- matrix(c(1, 1/sqrt(2), r/sqrt(2),
                   1/sqrt(2), 1, r,
                   r/sqrt(2), r, 1), nrow = 3)
  res <- as.numeric( -2*fMultivar::pnorm2d(-de1, 0, rho = 1/sqrt(2)) +
                       4*mnormt::pmnorm(c(-de1, 0, 0), mean = rep(0, 3), varcov = mat2) )
  return(res)
}
#' @rdname bridge
#' @import stats
#' @importFrom fMultivar pnorm2d
#' @importFrom mnormt pmnorm
bridgeF_ct <- function(r, zratio1 = NULL, zratio2){
  # continuous and truncated
  de1 <- qnorm(zratio2)
  mat2 <- matrix(c(1, 1/sqrt(2), r/sqrt(2),
                   1/sqrt(2), 1, r,
                   r/sqrt(2), r, 1), nrow = 3)
  res <- as.numeric( -2*fMultivar::pnorm2d(-de1, 0, rho = 1/sqrt(2)) +
                       4*mnormt::pmnorm(c(-de1, 0, 0), mean = rep(0, 3), varcov = mat2) )
  return(res)
}
#' @rdname bridge
#' @import stats
#' @importFrom fMultivar pnorm2d
#' @importFrom mnormt pmnorm
bridgeF_bb <- function(r, zratio1, zratio2){
  # binary and binary
  de1 <- qnorm(zratio1)
  de2 <- qnorm(zratio2)
  res <- as.numeric(2*(fMultivar::pnorm2d(de1, de2, rho = r) - zratio1*zratio2))
  return(res)
}
#' @rdname bridge
#' @import stats
#' @importFrom fMultivar pnorm2d
#' @importFrom mnormt pmnorm
bridgeF_tt <- function(r, zratio1, zratio2){
  # truncated and truncated
  de1 <- qnorm(zratio1)
  de2 <- qnorm(zratio2)

  mat1 <- matrix(c(1, 0, 1/sqrt(2), -r/sqrt(2),
                   0, 1, -r/sqrt(2), 1/sqrt(2),
                   1/sqrt(2), -r/sqrt(2), 1, -r,
                   -r/sqrt(2), 1/sqrt(2), -r, 1), nrow = 4)
  mat2 <- matrix(c(1, r, 1/sqrt(2), r/sqrt(2),
                   r, 1, r/sqrt(2), 1/sqrt(2),
                   1/sqrt(2), r/sqrt(2), 1, r,
                   r/sqrt(2), 1/sqrt(2), r, 1), nrow = 4)

  res <- as.numeric( -2*mnormt::pmnorm(c(-de1, -de2, 0, 0), mean = rep(0, 4), varcov = mat1) +
                       2*mnormt::pmnorm(c(-de1, -de2, 0, 0), mean = rep(0, 4), varcov = mat2)
  )
  return(res)
}
