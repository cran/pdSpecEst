## Basis Hermitian matrices

E_basis <- function(d) {
  .Deprecated("pdSpecEst:::E_coeff or pdSpecEst:::E_coeff_inv")
  E <- function(i, j) {
    E_ij <- matrix(0, nrow = d, ncol = d)
    if (i == j) {
      E_ij[i, j] <- 1
    } else if (i < j) {
      E_ij[c((j - 1) * d + i, (i - 1) * d + j)] <- 1/sqrt(2)
    } else {
      E_ij[c((j - 1) * d + i, (i - 1) * d + j)] <- complex(imaginary = c(-1, 1))/sqrt(2)
    }
    E_ij
  }

  indices <- expand.grid(1:d, 1:d)
  return(mapply(E, indices$Var1, indices$Var2, SIMPLIFY = "array"))
}

## Basis Tangent space

T_basis <- function(E, y) {
  .Deprecated("pdSpecEst:::T_coeff or pdSpecEst:::T_coeff_inv")
  d <- nrow(E)
  y.sqrt <- Sqrt(y)
  return(array(c(apply(E, 3, function(E) (y.sqrt %*% E) %*% y.sqrt)), dim = c(d, d, d^2)))
}



