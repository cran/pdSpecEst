
<!-- README.md is generated from README.Rmd. Please edit that file -->
The `pdSpecEst` package
=======================

[![CRAN version](http://www.r-pkg.org/badges/version/pdSpecEst)](https://cran.r-project.org/package=pdSpecEst) [![Travis-CI Build Status](https://travis-ci.org/JorisChau/pdSpecEst.svg?branch=master)](https://travis-ci.org/JorisChau/pdSpecEst)

The `pdSpecEst` (**p**ositive **d**efinite **Spec**tral **Est**imation) package provides data analysis tools for samples of symmetric or Hermitian positive definite matrices, such as collections of positive definite covariance matrices or spectral density matrices.

The tools in this package can be used to perform:

-   *Intrinsic wavelet transforms* for curves (1D) or surfaces (2D) of Hermitian positive definite matrices, with applications to dimension reduction, denoising and clustering for curves or surfaces of Hermitian positive definite matrices such as (time-varying) Fourier spectral density matrices. These implementations are based in part on the paper (Chau and Sachs 2017) and the material in Chapters 3 and 5 of (Chau 2018).

-   Exploratory data analysis and inference for samples of Hermitian positive definite matrices by means of *intrinsic data depth functions* and *depth rank-based hypothesis tests*. These implementations are based on the paper (Chau, Ombao, and Sachs 2017).

For more details and examples on how to use the package see the accompanying vignettes in the vignettes folder.

An R-Shiny app to demonstrate and test the implemented functionality in the package is available [here](https://jchau.shinyapps.io/pdSpecEst/).

*Author and maintainer:* Joris Chau (<j.chau@uclouvain.be>).

Installation
------------

-   **Stable CRAN version:** install from within R

-   **Current development version:** install via `devtools::install_github("JorisChau/pdSpecEst")`

References
----------

Chau, J. 2018. “Advances in Spectral Analysis for Multivariate, Nonstationary and Replicated Time Series.” PhD thesis, Université catholique de Louvain.

Chau, J., and R. von Sachs. 2017. “Intrinsic Wavelet Regression for Curves of Hermitian Positive Definite Matrices.” *ArXiv Preprint 1701.03314*. <https://arxiv.org/abs/1701.03314>.

Chau, J., H. Ombao, and R. von Sachs. 2017. “Intrinsic Data Depth for Hermitian Positive Definite Matrices.” *ArXiv Preprint 1706.08289*. <https://arxiv.org/abs/1706.08289>.
