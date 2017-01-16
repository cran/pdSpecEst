Impute_man <- function(M_j0, D, Nw) {
  d <- dim(M_j0)[1]
  j0 <- log(dim(M_j0)[3], 2)
  tM_j1 <- array(dim = c(d, d, 2^j0))
  P0 <- M_j0[, , 2^(j0 - 1)]
  M_j0.log <- sapply(1:2^j0, function(i) Logm(P0, M_j0[, , i]), simplify = "array")
  if (D == 0) {
    for (i in 1:2^(j0 - 1)) {
      P0 <- M_j0[, , 2 * i - 1]
      tM_j1[, , 2 * i - 1] <- Expm(P0, (1/4) * Logm(P0, M_j0[, , 2 * i]))
      tM_j1[, , 2 * i] <- Expm(P0, (5/4) * Logm(P0, M_j0[, , 2 * i]))
    }
  } else {
    for (i in 1:2^j0) {
      if (i == 1) {
        tM_j1[, , i] <- Expm(M_j0[, , i], (1/4) * Logm(M_j0[, , i], M_j0[, , i + 1]))
      } else if (i == 2^j0) {
        tM_j1[, , i] <- Expm(M_j0[, , i], (-1/4) * Logm(M_j0[, , i], M_j0[, , i - 1]))
      } else {
        nbrs0 <- abs((0:2^j0) - (0:2^j0)[i]) <= D
        Di <- min(sum(nbrs0[1:(i - 1)]), sum(nbrs0[(i + 1):2^j0]))
        nbrs <- abs((0:2^j0) - (0:2^j0)[i]) <= Di
        tM_j1[, , i] <- Expm(P0, apply(array(rep(Nw[[Di + 1]], each = d^2),
                              dim = c(d, d, 2 * Di + 1)) * M_j0.log[, , (i - Di):(i + Di)], c(1, 2), sum))
      }
    }
  }
  return(tM_j1)
}
