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



