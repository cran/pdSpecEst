# pdSpecEst 1.2.6

* Fixed CRAN check warnings/notes

# pdSpecEst 1.2.5

* Fixed CRAN check error

# pdSpecEst 1.2.4

* Fixed testthat error CRAN testing against R's messages
* Updated article references

# pdSpecEst 1.2.3

* Updated functions:
  + new metric `Riemannian-Rahman` in `WavTransf1D` and `InvWavTransf1D`
  + additional example spectra with `rExamples1D`

# pdSpecEst 1.2.3

* Fixed warning non-ASCII characters in references
* Fixed note CRAN checks

# pdSpecEst 1.2.2

* New and updated features:
  + intrinsic k-means clustering with `pdkMeans`
  + new example simulated data in `rExamples1D` and `rExamples2D`
  + intrinsic median of HPD matrices with `pdMedian`
* New and updated functions: `pdSpecClust1D`, `pdSpecClust2D`, `pdDepth`, `pdSpecEst1D`, `pdSpecEst2D`, `rExamples1D`, `rExamples2D`, `InvWavTransf1D`, `InvWavTransf2D`, `pdDist`, 
`pdMean`, `pdMedian`, `pdPgram`, `pdPgram2D`, `pdNeville`, `WavTransf1D`, `WavTransf2D`, 
`pdParTrans`, `Expm`, `Logm`
* Updated vignettes: 
  + "Wavelet-based multivariate Fourier spectral estimation", 
  + "Data depth and rank-based tests for HPD matrices"
* Updated Shiny app (see README or DESCRIPTION for url)

# pdSpecEst 1.2.1

* Removed unnecessary package imports

# pdSpecEst 1.2.0

* New and updated features: 
  + Updated 1D intrinsic wavelet regression and clustering with `pdSpecEst1D` and `pdSpecClust1D` based on various new metrics
  + 2D intrinsic wavelet regression and clustering with `pdSpecEst2D` and `pdSpecClust2D`
  + New tools for intrinsic 1D and 2D polynomial generation and interpolation
  + Time-varying periodograms with `pdPgram2D`
  + Tree-structured 1D and 2D wavelet thresholding with `pdCART` 
  + Depth-based intrinsic confidence regions with `pdConfInt1D` 
  + New example spectral matrices and benchmark procedures
* New and updated functions: `H.coeff`, `InvWavTransf1D`, `InvWavTransf2D`, `ParTrans`, `pdCART`, `pdConfInt1D`, `pdDepth`, `pdMean`, `pdNeville`, `pdPgram2D`, `pdPolynomial`, 
`pdRankTests`, `pdSpecClust1D`, `pdSpecClust2D`, `pdSpecEst1D`, `pdSpecEst2D`, `pdSplineReg`, `rExamples`, `rExamples2D`, `WavTransf1D`, `WavTransf2D`
* Updated vignettes: "Wavelet-based multivariate spectral analysis"
* Updated Shiny app (see README or DESCRIPTION for url)

# pdSpEcEst 1.1.1

* New features: data depth and rank-based hypothesis tests for samples of Hermitian PD matrices
* New functions: `pdDist`, `pdDepth`, `pdRankTests`
* New vignette: "Data depth and rank-based tests for HPD matrices"
* New demo Shiny app (see README or DESCRIPTION for url)

# pdSpecEst 1.0.0

* New release



