context("Spectral estimation and clustering")

test_that("Correctly working spectral estimation and clustering", {

  n <- 2^7
  d <- 2
  S <- 5
  ts.sim <- rARMA(n, d, array(0, c(d, d, 2)), array(0, c(d, d, 2)), diag(d))
  pgram <- pdPgram(ts.sim$X)
  f <- pdSpecEst(pgram$P)
  P_s <- replicate(S, pgram$P)
  cl <- pdSpecClust(P_s, K = 2, lam = f$lam)

  expect_equal(names(ts.sim), c("X", "f"))
  expect_equal(dim(ts.sim$X), c(n, d))
  expect_equal(InvWavTransf(WavTransf(pgram$P)$D), pgram$P)
  expect_equal(names(pgram), c("freq", "P"))
  expect_equal(dim(pgram$P), c(d, d, n/4))
  expect_equal(names(f), c("f", "D", "lam", "components"))
  expect_equal(dim(f$f), c(d, d, n/4))
  expect_equal(length(f$D), log2(n/4))
  expect_true(is.numeric(f$lam))
  expect_equal(length(f$components$not_thresholded), log2(n/4) - 1)
  expect_equal(length(f$components$thresholded), log2(n/4) - 1)
  expect_equal(dim(cl), c(S, 2))
  expect_equal(sum(cl), S)

})



