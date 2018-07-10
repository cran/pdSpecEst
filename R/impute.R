## Average-interpolation refinement scheme 2D
impute2D_R <- function(M_j0, L, metric, method){

  ## Set variables
  n <- dim(M_j0)[c(3,4)]
  d <- dim(M_j0)[1]
  L <- ifelse((2 * L + 1) > n, floor((n - 1) / 2), L)

  if(isTRUE(all.equal(L, c(0, 0)))){
    ## Trivial prediction
    tM1 <- M_j0[, , rep(1:n[1], each = 2), rep(1:n[2], each = 2)]

  } else {
    ## Nontrivial prediction using available weights
    M0 <- array(c(M_j0), dim = c(d, d, n[1] * n[2]))
    W <- lapply(1:length(W_2D), function(i) array(c(aperm(W_2D[[i]], c(3, 4, 1, 2))),
                                                dim = c(dim(W_2D[[i]])[3] * dim(W_2D[[i]])[4], 4)))
    tM1 <- impute2D_C(M0, W, n[1], n[2], L, metric, method)
    tM1 <- array(c(tM1), dim = c(d, d, 2 * n[1], 2 * n[2]))
  }
  return(tM1)
}
