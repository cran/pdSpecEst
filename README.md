
<!-- README.md is generated from README.Rmd. Please edit that file -->
The `pdSpecEst` package
=======================

[![CRAN version](http://www.r-pkg.org/badges/version/pdSpecEst)](https://cran.r-project.org/package=pdSpecEst) [![Travis-CI Build Status](https://travis-ci.org/JorisChau/pdSpecEst.svg?branch=master)](https://travis-ci.org/JorisChau/pdSpecEst)

The `pdSpecEst` (**p**ositive **d**efinite **Spec**tral **Est**imation) package provides data analysis tools for samples of symmetric or Hermitian positive definite matrices, such as collections of (non-degenerate) covariance matrices or spectral density matrices.

The tools in this package can be used to perform:

-   *Intrinsic manifold wavelet regression* and *clustering* for curves (1D) or surfaces (2D) of Hermitian positive definite matrices. These implementations are based in part on the paper (Chau and von Sachs 2017).

-   Exploratory data analysis and inference for samples of Hermitian positive definite matrices by means of *intrinsic manifold data depth* and *manifold rank-based hypothesis tests*. These implementations are based on the paper (Chau, Ombao, and von Sachs 2017).

For more details and examples on how to use the package see the accompanying vignettes in the vignettes folder.

A demo Shiny app to test out the implemented functions in the package is available [here](https://jchau.shinyapps.io/pdSpecEst/).

*Author and maintainer:* Joris Chau (<j.chau@uclouvain.be>).

Installation
------------

-   **Stable CRAN version:** install from within R

-   **Current development version:** install via `devtools::install_github("JorisChau/pdSpecEst")`

References
----------

Chau, J., and R. von Sachs. 2017. “Positive Definite Multivariate Spectral Estimation: A Geometric Wavelet Approach.” <http://arxiv.org/abs/1701.03314>.

Chau, J., H. Ombao, and R. von Sachs. 2017. “Data Depth and Rank-Based Tests for Covariance and Spectral Density Matrices.” <http://arxiv.org/abs/1706.08289>.
