<!-- README.md is generated from README.Rmd. Please edit that file -->
mixedCCA: sparse CCA for data of mixed types
============================================

The R package `mixedCCA` implements sparse canonical correlation analysis for data of mixed types: continuous, binary or zero-inflated (truncated continuous). The corresponding reference is

"Sparse semiparametric canonical correlation analysis for data of mixed types" by Yoon, Carroll and Gaynanova (2018+).

Installation
------------

``` install
library(devtools)
devtools::install_github("irinagain/mixedCCA")
```

Example
-------

``` r
library(mixedCCA)
#> Loading required package: MASS

# Data generation
set.seed(1)
n <- 100; p1 <- 5; p2 <- 7
maxcancor <- 0.9 # true canonical correlation
Sigma1 <- autocor(p1, 0.7)
groupind <- c(rep(1, 2), rep(2, p2-2))
Sigma2 <- blockcor(0.7, groupind)
copula1 <- "exp"; copula2 <- "cube" # transformation for Gaussian copula model
mu <- rep(0, p1+p2)
type1 <- type2 <- "trunc" # X1: truncated continuous, X2: truncated continuous
c1 <- rep(0, p1); c2 <- rep(0, p2) # threshold for truncation of underlying continuous variable
trueidx1 <- c(0, 0, 1, 1, 1) # true variable indices for dataset X1
trueidx2 <- c(1, 0, 1, 0, 0, 1, 1) # true variable indices for dataset X2
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2, maxcancor = maxcancor,
                       Sigma1 = Sigma1, Sigma2 = Sigma2,
                       copula1 = copula1, copula2 = copula2,
                       muZ = mu,
                       type1 = type1, type2 = type2, c1 = c1, c2 = c2
)
X1 <- simdata$X1
X2 <- simdata$X2

# estimated latent correlation matrix R
estimateR(X1, type = "trunc")
estimateR_mixed(X1, X2, type1 = "trunc", type2 = "trunc")


# sparse semiparametric CCA with BIC1 criterion
w1init <- rep(1, p1)
w2init <- rep(1, p2)
mixedCCAresult <- mixedCCA(X1, X2, type1 = "trunc", type2 = "trunc", lam.eps = 0.01, nlambda = 20,
                        w1init = w1init, w2init = w2init, BICtype = 1)
mixedCCAresult$KendallR # can extract latent correlation matrix estimated within the function
mixedCCAresult$w1 # canonical direction of X1
mixedCCAresult$w2 # canonical direction of X2
mixedCCAresult$cor # canonical correlation
```
