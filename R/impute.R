## Average-interpolation refinement scheme 1D
AIRefine1D <- function(M_j0, D, method, inverse = T, metric) {

  metric <- match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "rootEuclidean", "Euclidean"))
  n <- dim(M_j0)[3]
  d <- dim(M_j0)[1]
  tM_j1 <- (if(inverse) array(dim = c(d, d, 2 * n)) else array(dim = c(d, d, n)))
  D <- ifelse((2 * D + 1) > n, floor((n - 1) / 2), D)
  N <- 2 * D + 1

  if(missing(method)){
    method <- ifelse(N <= 9, "weights", "neville")
  }

  if(D == 0 | n == 1) {

    tM_j1 <- (if(inverse) M_j0[, , rep(1:n, each = 2)] else M_j0)

  } else {

    for(k in 1:n){

      ## Predict with available weights
      nbrs <- (if((k - D) < 1) 1:N else if((k + D) >= n) (n - N + 1):n else (k - D):(k + D))
      M <- M_j0[, , nbrs]

      if((N <= 9) & (method == "weights")){

        l <- (if((k - D) < 1) 2 * k else if((k + D) >= n) 2 * (N - (n - k)) else N + 1)
        w <- W_1D[[D + 1]][l, ]

        if(metric == "Riemannian"){

          Mean_new <- M_j0[, , k]
          Mean <- diag(d)
          m <- 0
          while((pdDist(Mean_new, Mean) > 1E-10) & (m < 20)){
            Mean <- Mean_new
            Mean_new <- Expm(Mean, apply(sapply(1:dim(M)[3], function(i) w[i] *
                                                  Logm(Mean, M[,,i]), simplify = "array"), c(1, 2), sum))
            m <- m+1
          }

        } else {
          Mean_new <- pdMean(M, w, metric = "Euclidean")
        }

        if(inverse){
          tM_j1[, , 2 * k] <- Mean_new
          tM_j1[, , 2 * k - 1] <- (if(metric == "Riemannian"){
            (M_j0[, , k] %*% solve(tM_j1[, , 2 * k])) %*% M_j0[, , k]
          } else 2 * M_j0[, , k] - tM_j1[, , 2 * k])
        } else{
          tM_j1[, , k] <- Mean_new
        }
      } else if(method == "neville"){

        M.bar <- sapply(1:N, function(i){
          pdMean(array(c(M[, , 1:i]), dim = c(d, d, i)),
                 metric = ifelse(metric == "Riemannian", "Riemannian", "Euclidean"))
        }, simplify = "array")

        if((k - D) < 1){
          ## At the left boundary
          tm1 <- (if(metric == "Riemannian"){
            Expm(M.bar[, , k], -(2 * k - 1) * Logm(M.bar[, , k],
                                                   pdNeville(P = M.bar, X = 1:N, x = k - 0.5)))
          } else {
            M.bar[, , k] - (2 * k - 1) * (pdNeville(P = M.bar, X = 1:N,
                                                    x = k - 0.5, metric = "Euclidean") - M.bar[, , k])
          })
        } else if ((k + D) >= n){
          ## At the right boundary
          k1 <- N - (n - k)
          tm1 <- (if(metric == "Riemannian"){
            Expm(M.bar[, , k1], -(2 * k1 - 1) * Logm(M.bar[, , k1],
                                                     pdNeville(P = M.bar, X = 1:N, x = k1 - 0.5)))
          } else {
            M.bar[, , k1] - (2 * k1 - 1) * (pdNeville(P = M.bar, X = 1:N,
                                                      x = k1 - 0.5, metric = "Euclidean") - M.bar[, , k1])
          })
        } else{
          ## Away from the boundary
          tm1 <- (if(metric == "Riemannian"){
            Expm(M.bar[, , D + 1], -N * Logm(M.bar[, , D + 1],
                                             pdNeville(P = M.bar, X = 1:N, x = D + 0.5)))
          } else {
            M.bar[, , D + 1] - N * (pdNeville(P = M.bar, X = 1:N,
                                              x = D + 0.5, metric = "Euclidean") - M.bar[, , D + 1])
          })
        }
        if(inverse){
          tM_j1[, , 2 * k] <- tm1
          tM_j1[, , 2 * k - 1] <- (if(metric == "Riemannian"){
            (M_j0[, , k] %*% solve(tm1)) %*% M_j0[, , k]
          } else 2 * M_j0[, , k] - tm1)
        } else{
          tM_j1[, , k] <- tm1
        }
      }
    }
  }
  return(tM_j1)
}

## Average-interpolation refinement scheme 2D
AIRefine2D <- function(M_j0, D, method, metric){

  ## Set variables
  metric = ifelse(metric == "Riemannian", "Riemannian", "Euclidean")
  j01 = log(dim(M_j0)[3], 2)
  j02 = log(dim(M_j0)[4], 2)
  d = dim(M_j0)[1]

  if(identical(j01, 0) & identical(j02, 0)){

    tM_j1 <- array(M_j0[, , 1, 1], dim = c(d, d, 2, 2))

  } else {

    tM_j1 <- array(dim = c(d, d, 2^(j01 + 1), 2^(j02 + 1)))
    grid_k <- expand.grid(1:2^j01, 1:2^j02)

    for(i in 1:nrow(grid_k)){

      D_k <- c(min(grid_k[i, 1] - 1, D[1]), min(grid_k[i, 2] - 1, D[2]))
      nbrs <- list(x = max(grid_k[i, 1] - D[1], 1):min(grid_k[i, 1] + D[1], 2^j01),
                   y = max(grid_k[i, 2] - D[2], 1):min(grid_k[i, 2] + D[2], 2^j02))
      N_k <- 2 * D_k + 1

      if((length(nbrs$x) == 2 * D[1] + 1) & (length(nbrs$y) == 2 * D[2] + 1) & all(D <= 4) & (method == "weights")){

        ## Away from the boundary:
        ## Predict with available weights
        pred2D <- function(l1, l2){
          M <- array(c(M_j0[, , nbrs$x, nbrs$y]), dim = c(d, d, length(nbrs$x) * length(nbrs$y)))
          w <- c(W_2D[[5 * D_k[2] + (D_k[1] + 1)]][l1, l2, , ])

          if(metric == "Riemannian"){
            Mean_new <- M_j0[, , grid_k[i, 1], grid_k[i, 2]]
            Mean <- diag(d)
            m <- 0
            while((pdDist(Mean_new, Mean) > 1E-8) & (m < 15)){
              Mean <- Mean_new
              Mean_new <- Expm(Mean, apply(sapply(1:dim(M)[3], function(k) w[k] *
                                                    Logm(Mean, M[, , k]), simplify = "array"), c(1, 2), sum))
              m <- m+1
            }
          } else{
            Mean_new <- pdMean(M, w, metric = "Euclidean")
          }
          return(Mean_new)
        }

        tM_j1[, , (2 * grid_k[i, 1] - 1):(2 * grid_k[i, 1]), (2 * grid_k[i, 2] - 1):(2 * grid_k[i, 2])] <-
          array(c(mapply(pred2D, rep(c(1, 2), times = 2), rep(c(1, 2), each = 2), SIMPLIFY = "array")), dim = c(d, d, 2, 2))

      } else {

        ## At the boundary or if max(order) > 9,
        ## Predict with Neville's algorithm
        grid_l <- expand.grid(nbrs$x, nbrs$y)

        M.bar <- array(c(mapply(function(i1,i2){
          pdMean(array(c(M_j0[, , grid_l[1, 1]:i1, grid_l[1, 2]:i2]),
                       dim = c(d, d, length(grid_l[1, 1]:i1) * length(grid_l[1, 2]:i2))), metric = metric)
        }, grid_l$Var1, grid_l$Var2, SIMPLIFY="array")), dim = c(d, d, length(nbrs$x), length(nbrs$y)))

        M.bar_pred <- pdNeville(P = M.bar, X = nbrs, x = list(x = grid_k[i, 1] - 1 + c(0, 0.5, 1),
                                                              y = grid_k[i, 2] - 1 + c(0, 0.5, 1)), metric = metric)

        ## Predict M_{j+1,2k_1+1,2k_2+1}
        C0 <- pdMean(array(c(M.bar_pred[, , 3, 2], M.bar_pred[, , 2, 2]), dim=c(d, d, 2)),
                     w = c(N_k[1] + 1, -N_k[1]), metric = metric)
        M.bar_min0 <- pdMean(array(c(M.bar_pred[, , 2, 3], C0), dim=c(d, d, 2)), metric = metric,
                             w = c((N_k[1] * (N_k[2] + 1)) / ((N_k[1] + 1) * (N_k[2] + 1) - 1), N_k[2] / ((N_k[1] + 1) * (N_k[2] + 1) - 1)))
        tM_j1[, , 2 * grid_k[i, 1], 2 * grid_k[i, 2]] <- pdMean(array(c(M.bar[, , which(nbrs$x == grid_k[i, 1]),
                                                                which(nbrs$y == grid_k[i, 2])], M.bar_min0), dim = c(d, d, 2)), metric = metric,
                                                                w = c((N_k[1] + 1) * N_k[2] + N_k[1] + 1, -((N_k[1] + 1) * N_k[2] + N_k[1])))

        ## Predict M_{j+1,2k_1,2k_2+1}
        C1 <- pdMean(array(c(M.bar_pred[, , 1, 3], M.bar_pred[, , 1, 2]), dim = c(d, d, 2)), w = c(N_k[2] + 1, -N_k[2]), metric = metric)
        M.bar_min1 <- pdMean(array(c(M.bar_pred[, , 3, 2], C1), dim=c(d, d, 2)), metric = metric,
                             w = c((N_k[1] + 1) * N_k[2] / ((N_k[1] + 1) * (N_k[2] + 1) - 2), (N_k[1] - 1) / ((N_k[1] + 1) * (N_k[2] + 1) - 2)))
        M.block1 <- pdMean(array(c(M.bar[, , which(nbrs$x == grid_k[i, 1]), which(nbrs$y == grid_k[i, 2])], M.bar_min1), dim = c(d, d, 2)),
                           w = c(((N_k[1] + 1) * (N_k[2] + 1) - 2) / 2 + 1, -((N_k[1] + 1) * (N_k[2] + 1) - 2) / 2), metric = metric)
        tM_j1[, , 2 * grid_k[i, 1] - 1, 2 * grid_k[i, 2]] <- pdMean(array(c(M.block1, tM_j1[, , 2 * grid_k[i, 1], 2 * grid_k[i, 2]]),
                                                                          dim = c(d, d, 2)), w = c(2, -1), metric = metric)

        ## Predict M_{j+1,2k_1+1,2k_2}
        C2 <- pdMean(array(c(M.bar_pred[, , 3, 1], M.bar_pred[, , 2, 1]), dim = c(d, d, 2)), w = c(N_k[1] + 1, -N_k[1]), metric = metric)
        M.bar_min2 <- pdMean(array(c(M.bar_pred[, , 2, 3], C2), dim = c(d, d, 2)), metric = metric,
                             w = c((N_k[2] + 1) * N_k[1] / ((N_k[2] + 1) * (N_k[1] + 1) - 2), (N_k[2] - 1) / ((N_k[2] + 1) * (N_k[1] + 1) - 2)))
        M.block2 <- pdMean(array(c(M.bar[, , which(nbrs$x == grid_k[i, 1]), which(nbrs$y == grid_k[i, 2])], M.bar_min2), dim = c(d, d, 2)),
                           w = c(((N_k[1] + 1) * (N_k[2] + 1) - 2) / 2 + 1, -((N_k[1] + 1) * (N_k[2] + 1) - 2) / 2), metric = metric)
        tM_j1[, , 2 * grid_k[i, 1], 2 * grid_k[i, 2] - 1] <- pdMean(array(c(M.block2, tM_j1[, , 2 * grid_k[i, 1], 2 * grid_k[i, 2]]),
                                                                          dim = c(d, d, 2)), w = c(2, -1), metric = metric)

        ## Predict M_{j+1,2k_1,2k_2}
        if(metric == "Riemannian"){
          tM_j1[, , 2 * grid_k[i, 1] - 1, 2 * grid_k[i, 2] - 1] <- Expm(M_j0[, , grid_k[i, 1], grid_k[i, 2]],
                                                                        -(Logm(M_j0[, , grid_k[i, 1], grid_k[i, 2]], tM_j1[, , 2 * grid_k[i, 1], 2 * grid_k[i, 2] - 1]) +
                                                                            Logm(M_j0[, , grid_k[i, 1], grid_k[i, 2]], tM_j1[, , 2 * grid_k[i, 1] - 1, 2 * grid_k[i, 2]]) +
                                                                            Logm(M_j0[, , grid_k[i, 1], grid_k[i, 2]], tM_j1[, , 2 * grid_k[i, 1], 2 * grid_k[i, 2]])))
        } else{
          tM_j1[, , 2 * grid_k[i, 1] - 1, 2 * grid_k[i, 2] - 1] <- 4 * M_j0[, , grid_k[i, 1], grid_k[i, 2]] -
                                                                            tM_j1[, , 2 * grid_k[i, 1], 2 * grid_k[i, 2] - 1] -
                                                                            tM_j1[, , 2 * grid_k[i, 1] - 1, 2 * grid_k[i, 2]] -
                                                                            tM_j1[, , 2 * grid_k[i, 1], 2 * grid_k[i, 2]]
        }
      }
    }
  }

  return(tM_j1)
}

