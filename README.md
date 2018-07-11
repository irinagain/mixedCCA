<!-- README.md is generated from README.Rmd. Please edit that file -->
mixedCCA: sparse CCA for data of mixed types
============================================

The R package `mixedCCA` is to apply sparse CCA on the data of mixed types: continuous, binary or zero-inflated (truncated continuous). The methods are described in "Sparse semiparametric canonical correlation analysis for data of mixed types" by Yoon, Carroll and Gaynanova (2018+).

Installation
------------

``` install
library(devtools)
devtools::install_github("irinagain/mixedCCA")
```

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
## basic example code
library(mixedCCA)
#> Loading required package: MASS

# Data generation
set.seed(1)
n <- 100; p1 <- 5; p2 <- 7
maxcancor <- 0.9 # true canonical correlation
Sigma1 <- autocor(p1, 0.7)
groupind <- c(rep(1, 2), rep(2, p2-2))
Sigma2 <- blockcor(0.7, groupind)
copula1 <- "exp"; copula2 <- "cube" # transformation
mu <- rep(0, p1+p2)
type1 <- type2 <- "trunc" # X1: truncated continuous, X2: truncated continuous
c1 <- rep(0, p1); c2 <- rep(0, p2)
trueidx1 <- c(0, 0, 1, 1, 1)
trueidx2 <- c(1, 0, 1, 0, 0, 1, 1)
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2, maxcancor = maxcancor,
                       Sigma1 = Sigma1, Sigma2 = Sigma2,
                       copula1 = copula1, copula2 = copula2,
                       muZ = mu,
                       type1 = type1, type2 = type2, c1 = c1, c2 = c2
)

X1 <- simdata$X1
X2 <- simdata$X2

zratio1 <- colMeans(X1==0)
zratio2 <- colMeans(X2==0)

bridge <- bridge_select(type1 = "trunc", type2 = "trunc")
K12 <- Kendall_matrix(X1, X2)

# estimated latent correlation R based on K12.
fromKtoR_mixed(K12, zratio1 = zratio1, zratio2 = zratio2, type1 = "trunc", type2 = "trunc")

# 
estimateR(X1, type = "trunc")
estimateR_mixed(X1, X2, type1 = "trunc", type2 = "trunc")
```
