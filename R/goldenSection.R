gss <- function(range, cv, tol = 0.01) {
  x1 <- range[1]
  x4 <- range[2]
  x2 <- x4 + (x1 - x4)/((1 + sqrt(5))/2)
  x3 <- x1 + (x4 - x1)/((1 + sqrt(5))/2)
  fx2 <- cv(x2)
  fx3 <- cv(x3)
  i <- 0
  while (!isTRUE(all.equal(fx2, fx3)) & (abs(x2 - x3) > tol) & (i <= 20)) {
    i <- i + 1
    if (fx2 < fx3) {
      x4 <- x3
      x3 <- x2
      fx3 <- fx2
      x2 <- x4 + (x1 - x4)/((1 + sqrt(5))/2)
      fx2 <- cv(x2)
    } else {
      x1 <- x2
      x2 <- x3
      fx2 <- fx3
      x3 <- x1 + (x4 - x1)/((1 + sqrt(5))/2)
      fx3 <- cv(x3)
    }
  }
  return(mean(c(x2, x3)))
}
