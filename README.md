
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BDcocolasso

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/celiaescribe/BDcocolasso.svg?branch=master)](https://travis-ci.org/celiaescribe/BDcocolasso)
<!-- badges: end -->

R software package to implement high-dimensional error-in-variables
regression. This package implements CoCoLasso algorithm in settings with
additive error or missing data in the covariates. This package also
implements a variation of the CoCoLasso algorithm called Block-Descent
CoCoLasso (or BD-CoCoLasso), which focuses on a setting where only a
small percentage of the features are corrupted (with additive error or
missing data)

This package is based on the [CoCoLasso
algorithm](https://arxiv.org/pdf/1510.07123.pdf). CoCoLASSO requires a
computationally demanding positive semi-definite projection of the
covariance matrix for a high dimensional feature set. In a very
high-dimensional context where there are both corrupted and uncorrupted
covariates and where the portion of corrupted features is relatively small,
we take advantage of the block descent minimization to develop a
more efficient algorithm called BDCoCoLasso. In an alternating block
minimization algorithm, the CoCoLasso corrections are used when updating
corrupted coefficient vectors, and a simple LASSO is used for the
uncorrupted coefficient vectors. Both sub-problems are convex and hence a
global solution can be obtained, even though adaption of the
cross-validation step requires care in this setting where there are
products of corrupted and uncorrupted matrices.

## Installation

``` r
install.packages("BDcocolasso")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("celiaescribe/BDcocolasso")
```

## Vignette

See the online vignette for details about the BDcoco model and example
usage of the functions.

## Model input

There exist two settings in which the BD-CoCoLasso can be used : in the
simple CoCoLasso version, and in the Block-Descent-CoCoLasso version.
The inputs vary according to the chosen algorithm setting, and according
to the chosen noise setting.

  - **CoCoLasso setting**: This method requires seven inputs (let n be
    the number of observations and p the number of X variables):

<!-- end list -->

1.  **X**: n x p matrix of covariates. Can be high-dimensional, i.e., p
    \>\> n. Must be continuous or with binary categorical features. Can
    contain missing values in NA format in the missing data setting.
2.  **y**: a continuous response of length n.
3.  **n**: Number of samples
4.  **p**: Number of covariates
5.  **noise**: Type of noise setting. There are two possibilities :
    *additive* or *missing*. In the *additive* setting it is necessary
    to specify the *tau* parameter, corresponding to the standard
    deviation of the additive error matrix. In the *missing* setting,
    nothing has to be specified.
6.  **block**: Chosen setting. Here, *block* should be equal to *FALSE*.
7.  **penalty**: Type of penalty chosen. It can be equal to *lasso* or
    *SCAD* according to the chosen penalty setting.

<!-- end list -->

  - **BD-CoCoLasso setting**: This method requires nine inputs (let n be
    the number of observations, p the number of X variables, p1 the
    number of uncorrupted variables and p2 the number of corrupted
    variables, with p1 + p2 = p):

<!-- end list -->

1.  **X**: n x p matrix of covariates. Can be high-dimensional, i.e., p
    \>\> n. Must be continuous or with binary categorical features. Can
    contain missing values in NA format in the missing data setting. The
    first p1 columns must correspond to the uncorrupted covariates, and
    the last p2 columns must correspond to the corrupted covariates.
2.  **y**: a continuous response of length n.
3.  **n**: Number of samples
4.  **p**: Number of covariates
5.  **p1**: Number of uncorrupted covariates
6.  **p2**: Number of corrupted covariates
7.  **noise**: Type of noise setting. There are two possibilities :
    *additive* or *missing*. In the *additive* setting it is necessary
    to specify the *tau* parameter, corresponding to the standard
    deviation of the additive error matrix. In the *missing* setting,
    nothing has to be specified.
8.  **block**: Chosen setting. Here, *block* should be equal to *TRUE*.
9.  **penalty**: Type of penalty chosen. It can be equal to *lasso* or
    *SCAD* according to the chosen penalty setting.

<!-- end list -->

  - **Three-block BD-CoCoLasso setting**: This method handles a mixed error
    setting where both additive error and missing data occur. This requires 
    excecuting the function *generalcoco*. The required inputs are the same as
    the **BD-CoCoLasso** setting except that **p2** stands for the number of
    corrupted covariates measured with additive error and an additional parameter
    **p3** stands for the number of corrupted covariates measured with missing 
    data. It is essential in both settings that the covariates are sorted with
    the uncorrupted covariates are in the first columns. In the Three-block
    setting, the additive-error-containing covariates should precede the 
    missing-data-containing covariates as well.

## Contact

email : celia.escribe@polytechnique.edu

## Credit

We based this R package on the following articles :

  - [Cocolasso for high dimensional error-in-variables
    regression](https://arxiv.org/pdf/1510.07123.pdf)
  - [High-dimensional regression with noisy and missing data: Provable
    guarantees with nonconvexity](https://arxiv.org/pdf/1109.3714.pdf)
