<!-- README.md is generated from README.Rmd. Please edit that file -->
mixedCCA: sparse CCA for data of mixed types
============================================

The R package `mixedCCA` is to apply sparse CCA on the data of mixed types: continuous, binary or zero-inflated (truncated continuous). The methods are described in \[Sparse semiparametric canonical correlation analysis for data of mixed types\] by Yoon, Carroll and Gaynanova (2018+).

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
```
