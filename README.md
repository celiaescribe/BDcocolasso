
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BDcocolasso

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/celiaescribe/BDcocolasso.svg?branch=master)](https://travis-ci.org/celiaescribe/BDcocolasso)
<!-- badges: end -->

This package has two aims : the first one is to offer an implementation
of the CoCoLasso algorithm, and the second is to propose an implemtation
of the block descent cocolasso (BDCocolasso). CoCoLASSO requires a
computationally demanding positive semi-definite projection of the
covariance matrix for a high dimensional feature set. In our context
when there are corrupted and uncorrupted covariates, we take advantage
of the block descent minimization trick to develop a more efficient
algorithm. In an alternating block minimization algorithm, the CoCoLasso
corrections are used when updating corrupted coefficient vectors, and a
simple LASSO is used for the uncorrupted coefficient vectors. Both
subproblems are convex and hence a global solution can be obtained, even
though adaption of the cross-validation step requires care in this
setting where there are products of corrupted and uncorrupted matrices.

## Installation

``` r
install.packages("BDcocolasso")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("celiaescribe/BDcocolasso")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(BDcocolasso)
## basic example code
```
