context("Geometric functions")

test_that("Correctly working geometric tools", {

  d <- 2
  p1 <- matrix(complex(real = rnorm(d^2), imaginary = rnorm(d^2)), nrow = d)
  p2 <- matrix(complex(real = rnorm(d^2), imaginary = rnorm(d^2)), nrow = d)
  P1 <- t(Conj(p1)) %*% p1
  P2 <- t(Conj(p2)) %*% p2

  expect_equal(Expm(P1, Logm(P1, P2)), P2)
  expect_equal(c(Sqrt(P1) %*% iSqrt(P1)), as.complex(diag(d)))
  expect_equal(E_coeff_inv(E_coeff(P1)), P1)
  expect_equal(T_coeff_inv(T_coeff(P1, P2), P2), P1)
  expect_equal(KarchMean(array(c(P1, P2), dim = c(d, d, 2))), Mid(P1, P2))
  expect_true(is.numeric(pdDist(P1, P2)))
  expect_true(is.numeric(pdDist(P1, P2, method = "logEuclidean")))
  expect_true(is.numeric(pdDist(P1, P2, method = "Cholesky")))
  expect_true(is.numeric(pdDist(P1, P2, method = "Euclidean")))
  expect_true(is.numeric(pdDist(P1, P2, method = "Procrustes")))

})
