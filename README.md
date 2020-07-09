<!-- README.md is generated from README.Rmd. Please edit that file -->

mixedCCA: sparse CCA for data of mixed types
============================================

The R package `mixedCCA` implements sparse canonical correlation
analysis for data of mixed types: continuous, binary or zero-inflated
(truncated continuous). The corresponding reference is

[Yoon G., Carroll R.J. and Gaynanova I. (2020). “Sparse semiparametric
canonical correlation analysis for data of mixed types”.
*Biometrika*](https://doi.org/10.1093/biomet/asaa007).

The faster version of latent correlation computation part is now fully
available and implemented to the R package `mixedCCA`. The corresponding
reference is available on arXiv:

[Yoon G., Müller C.L. and Gaynanova I., “Fast computation of latent
correlations” *arXiv*](https://arxiv.org/abs/2006.13875).

Installation
------------

``` install
devtools::install_github("irinagain/mixedCCA")
```

Example
-------

``` r
library(mixedCCA)

### Simple example

# Data setting
n <- 100; p1 <- 15; p2 <- 10 # sample size and dimensions for two datasets.
maxcancor <- 0.9 # true canonical correlation

# Correlation structure within each data set
set.seed(0)
perm1 <- sample(1:p1, size = p1);
Sigma1 <- autocor(p1, 0.7)[perm1, perm1]
blockind <- sample(1:3, size = p2, replace = TRUE);
Sigma2 <- blockcor(blockind, 0.7)
mu <- rbinom(p1+p2, 1, 0.5)

# true variable indices for each dataset
trueidx1 <- c(rep(1, 3), rep(0, p1-3))
trueidx2 <- c(rep(1, 2), rep(0, p2-2))

# Data generation
simdata <- GenerateData(n=n, trueidx1 = trueidx1, trueidx2 = trueidx2, maxcancor = maxcancor,
                        Sigma1 = Sigma1, Sigma2 = Sigma2,
                        copula1 = "exp", copula2 = "cube",
                        muZ = mu,
                        type1 = "trunc", type2 = "trunc",
                        c1 = rep(1, p1), c2 =  rep(0, p2)
)
X1 <- simdata$X1
X2 <- simdata$X2

# Sparse semiparametric CCA with BIC1 criterion
mixedCCAresult <- mixedCCA(X1, X2, type1 = "trunc", type2 = "trunc", BICtype = 1)
mixedCCAresult$KendallR # estimated latent correlation matrix
mixedCCAresult$w1 # estimated canonical vector for X1
mixedCCAresult$w2 # estimated canonical vector for X2
mixedCCAresult$cancor # estimated canonical correlation

# Separate estimation of latent correlation matrix
estimateR(X1, type = "trunc")$R # For X1 only
estimateR_mixed(X1, X2, type1 = "trunc", type2 = "trunc")$R12 # For X = (X1, X2)
```
