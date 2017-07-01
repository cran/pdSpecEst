## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ------------------------------------------------------------------------
library(pdSpecEst)
## Fix parameters
freq <- seq(from = pi / 2^9, to = pi, length = 2^9)
d <- 2
Phi <- array(c(0.5, 0, 0, 0.2, 0, 0, 0, -0.9), dim = c(d, d, 2))
Theta <- array(c(0, 0.1, 0.1, 0, 0, 0, 0, 0.5), dim = c(d, d, 2))
Sigma <- matrix(c(2, 0, 0, 0.25), nrow = d)

## Generate time series
set.seed(0)
ts.sim <- rARMA(2^11, d, Phi, Theta, Sigma, freq = freq)
str(ts.sim)

## ---- eval=FALSE---------------------------------------------------------
#  ## Plot time series observations
#  par(mfrow = c(2,1), mar = c(4.5, 3, 2, 2))
#  invisible(sapply(1:d, function(i) plot(ts.sim$X[, i], main = paste0("Component ", i), type = "l", xlab = "Time", ylab = "")))
#  
#  ## Plot spectral matrix
#  layout(mat = matrix(c(1,1,2,3,4,5,6,6), nrow = 4))
#  plotspec <- function(i){
#    if(i[1] == i[2]){
#      plot(freq, Re(ts.sim$f[i[1], i[1], ]), main = paste0("Auto-spectrum (", i[1], ", ", i[1], ")"), type = "l", xlab = "Frequency", ylab = "")
#    } else{
#      plot(freq, Re(ts.sim$f[i[1], i[2], ]), main = paste0("Real cross-spectrum (", i[1], ", ", i[2], ")"), type = "l", xlab = "Frequency", ylab = "")
#      plot(freq, Im(ts.sim$f[i[1], i[2], ]), main = paste0("Imag. cross-spectrum (", i[1], ", ", i[2], ")"), type = "l", xlab = "Frequency", ylab = "")
#    }
#  }
#  invisible(apply(expand.grid(1:d, 1:d), 1, plotspec))

## ------------------------------------------------------------------------
pgram <- pdPgram(ts.sim$X)
str(pgram)

## ------------------------------------------------------------------------
f.hat <- pdSpecEst(pgram$P)
str(f.hat)

## ---- echo=FALSE, fig.show='hold', fig.cap='Figure: True (left) and estimated (right) spectral matrices in the frequency domain.', out.width='35%'----
## Plot spectral matrix
par(mar = c(0.1,0.1,2,0.1))
layout(mat = matrix(c(1,1,2,3,4,5,6,6), nrow = 4))
plotspec <- function(i, data, col){
  ylim <- range(Re(data$f), Im(data$f))
  if(i[1] == i[2]){
    plot(freq, Re(data$f[i[1], i[1], ]), main = paste0("Re(f(", i[1], ",", i[1], "))"), type = "l", xaxt = "n", yaxt = "n", col = col, ylim = ylim)
    abline(h = 0, lty = 3)
  } else{
    plot(freq, Re(data$f[i[1], i[2], ]), main = paste0("Re(f(", i[1], ",", i[2], "))"), type = "l", xaxt = "n", yaxt = "n", col = col, ylim = ylim)
    abline(h = 0, lty = 3)
    plot(freq, Im(data$f[i[1], i[2], ]), main = paste0("Im(f(", i[1], ",", i[1], "))"), type = "l", xaxt = "n", yaxt = "n", col = col, ylim = ylim)
    abline(h = 0, lty = 3)
  }
}
grid <- expand.grid(1:d, 1:d)
invisible(apply(grid, 1, function(i) plotspec(i, ts.sim, 1)))
invisible(apply(grid, 1, function(i) plotspec(i, f.hat, 2)))

## ------------------------------------------------------------------------
Phi1 <- array(c(0.5, 0, 0, 0.1, 0, 0, 0, -0.9), dim = c(d, d, 2))
Phi2 <- array(c(0.5, 0, 0, 0.3, 0, 0, 0, -0.9), dim = c(d, d, 2))
pgram <- function(Phi) pdPgram(rARMA(2^10, d, Phi, Theta, Sigma)$X)$P
P <- array(c(replicate(5, pgram(Phi1)), replicate(5, pgram(Phi2))), dim=c(d, d, 2^8, 10))
pdSpecClust(P, K = 2, lam = 3)

