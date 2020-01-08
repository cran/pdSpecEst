
<!-- README.md is generated from README.Rmd. Please edit that file -->

# The `pdSpecEst` package

The `pdSpecEst` (**p**ositive **d**efinite **Spec**tral **Est**imation)
package provides data analysis tools for samples of symmetric or
Hermitian positive definite matrices, such as collections of positive
definite covariance matrices or spectral density matrices.

The tools in this package can be used to perform:

  - *Intrinsic wavelet transforms* for curves (1D) or surfaces (2D) of
    Hermitian positive definite matrices, with applications to for
    instance: dimension reduction, denoising and clustering for curves
    or surfaces of Hermitian positive definite matrices such as
    (time-varying) Fourier spectral density matrices. These
    implementations are based in part on the papers (Chau and Sachs
    2019) and (Chau and Sachs 2018) and Chapters 3 and 5 of (Chau 2018).

  - Exploratory data analysis and inference for samples of Hermitian
    positive definite matrices by means of *intrinsic data depth
    functions* and *depth rank-based hypothesis tests*. These
    implementations are based on the paper (Chau, Ombao, and Sachs 2019)
    and Chapter 4 of (Chau 2018).

For more details and examples on how to use the package see the
accompanying vignettes in the vignettes folder.

*Author and maintainer:* Joris Chau (<joris.chau@openanalytics.eu>).

## Installation

  - **Stable CRAN version:** install from within R

## References

<div id="refs" class="references">

<div id="ref-C18">

Chau, J. 2018. “Advances in Spectral Analysis for Multivariate,
Nonstationary and Replicated Time Series.” PhD thesis, Universite
catholique de Louvain.

</div>

<div id="ref-COvS17">

Chau, J., H. Ombao, and R. von Sachs. 2019. “Intrinsic Data Depth for
Hermitian Positive Definite Matrices.” *Journal of Computational and
Graphical Statistics* 28 (2): 427–39.
<https://doi.org/https://doi.org/10.1080/10618600.2018.1537926>.

</div>

<div id="ref-CvS18">

Chau, J., and R. von Sachs. 2018. “Intrinsic Wavelet Regression for
Surfaces of Hermitian Positive Definite Matrices.” *ArXiv Preprint
1808.08764*. <https://arxiv.org/abs/1808.08764>.

</div>

<div id="ref-CvS17">

———. 2019. “Intrinsic Wavelet Regression for Curves of Hermitian
Positive Definite Matrices.” *Journal of the American Statistical
Association*.
<https://doi.org/https://doi.org/10.1080/01621459.2019.1700129>.

</div>

</div>
