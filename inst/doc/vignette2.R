## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 10, fig.height = 4, fig.align = "center", out.width = "600px")

## ------------------------------------------------------------------------
library(pdSpecEst)
set.seed(100)

## Pointwise random sample
X1 <- replicate(50, Expm(diag(2), pdSpecEst:::E_coeff_inv(0.5 * rnorm(4)))) 
str(X1)

## Curve random sample
X2 <- replicate(50, sapply(1:5, function(i) Expm(i * diag(2), pdSpecEst:::T_coeff_inv(0.5 * rnorm(4), i * diag(2))), simplify = "array"))
str(X2)

## ------------------------------------------------------------------------
## Pointwise depth
pdDepth(y = diag(2), X = X1, method = "gdd") ## geodesic distance depth
pdDepth(y = diag(2), X = X1, method = "zonoid") ## manifold zonoid depth
pdDepth(y = diag(2), X = X1, method = "spatial") ## manifold spatial depth

## Integrated depth 
pdDepth(y = sapply(1:5, function(i) i * diag(2), simplify = "array"), X = X2, method = "gdd") 
pdDepth(y = sapply(1:5, function(i) i * diag(2), simplify = "array"), X = X2, method = "zonoid") 
pdDepth(y = sapply(1:5, function(i) i * diag(2), simplify = "array"), X = X2, method = "spatial") 

## ------------------------------------------------------------------------
(dd1 <- pdDepth(X = X1, method = "gdd")) ## pointwise geodesic distance depth

(dd2 <- pdDepth(X = X2, method = "gdd")) ## integrated geodesic distance depth

## ------------------------------------------------------------------------
(dd1.ranks <- rank(1 - dd1)) ## pointwise depth ranks

(dd2.ranks <- rank(1 - dd2)) ## integrated depth ranks

## Explore sample X1
head(order(dd1.ranks)) ## most central observations 
rev(tail(order(dd1.ranks))) ## most outlying observations
X1[ , , which(dd1.ranks == 1)] ## most central HPD matrix 
X1[ , , which(dd1.ranks == 50)] ## most outlying HPD matrix


## ------------------------------------------------------------------------
(mean.X1 <- KarchMean(X1)) 

pdDepth(y = mean.X1, X = X1, method = "gdd")

## ---- echo = FALSE-------------------------------------------------------
d <- round(exp(seq(from = log(2), to = log(64), length = 10)), digits = 0)
n <- round(exp(seq(from = log(100), to = log(2000), length = 10)), digits = 0)
par(mfrow = c(1,2), mgp = c(2,0.5,0), mar = c(4,3,2,2))
plot(range(d), range(pdSpecEst:::depth_time, na.rm=T), type="n", log = "xy", main = "Increasing dimension d x d, (n = 500)", xlab = "d", ylab = "Time (ms)")
for(i in 1:3){
  lines(d, pdSpecEst:::depth_time[i,,1], col = i, lty = i)
  points(d, pdSpecEst:::depth_time[i,,1], col = i, pch = 15+i)
}
legend("bottomright", legend = c("gdd", "zonoid", "spatial"), lty = 1:3, pch = 16:18,
       col = 1:3, cex = 1, bty = "n", y.intersp = 1)
plot(range(n), range(pdSpecEst:::depth_time, na.rm=T), type="n", log = "xy", main = "Increasing sample size n, (d = 6)", xlab = "n", ylab = "Time (ms)")
for(i in 1:3){
  lines(n, pdSpecEst:::depth_time[i,,2], col = i, lty = i)
  points(n, pdSpecEst:::depth_time[i,,2], col = i, pch = 15+i)
}
legend("bottomright", legend = c("gdd", "zonoid", "spatial"), lty = 1:3, pch = 16:18,
       col = 1:3, cex = 1, bty = "n", y.intersp = 1)

## ------------------------------------------------------------------------
## Generate data (null true)
data1 <- array(c(X1, replicate(50, Expm(diag(2), pdSpecEst:::E_coeff_inv(0.5 * rnorm(4))))), dim = c(2, 2, 100)) ## pointwise sample
data2 <- array(c(X2, replicate(50, sapply(1:5, function(i) Expm(i * diag(2), pdSpecEst:::T_coeff_inv(0.5 * rnorm(4), i * diag(2))), simplify = "array"))), dim = c(2, 2, 5, 100)) ## curve sample

## Generate data (null false)
data1a <- array(c(X1, replicate(50, Expm(diag(2), pdSpecEst:::E_coeff_inv(rnorm(4))))), dim = c(2, 2, 100)) ## pointwise scale change
data2a <- array(c(X2, replicate(50, sapply(1:5, function(i) Expm(i * diag(2), pdSpecEst:::T_coeff_inv(rnorm(4), i * diag(2))), simplify = "arra"))), dim = c(2, 2, 5, 100)) ## curve scale change

## Rank-sum test
pdRankTests(data1, sample.sizes = c(50, 50), "rank.sum")[1:4] ## null true (pointwise)
pdRankTests(data2, sample.sizes = c(50, 50), "rank.sum")[2] ## null true (curve)
pdRankTests(data1a, sample.sizes = c(50, 50), "rank.sum")[2] ## null false (pointwise)
pdRankTests(data2a, sample.sizes = c(50, 50), "rank.sum")[2] ## null false (curve)

## Kruskal-Wallis test
pdRankTests(data1, sample.sizes = c(50, 25, 25), "krusk.wall")[1:4] ## null true (pointwise)
pdRankTests(data2, sample.sizes = c(50, 25, 25), "krusk.wall")[2] ## null true (curve)
pdRankTests(data1a, sample.sizes = c(50, 25, 25), "krusk.wall")[2] ## null false (pointwise)
pdRankTests(data2a, sample.sizes = c(50, 25, 25), "krusk.wall")[2] ## null false (curve)

## ------------------------------------------------------------------------
## Trial-specific means
mu <- replicate(50, Expm(diag(2), pdSpecEst:::E_coeff_inv(0.1 * rnorm(4))))

## Generate paired samples X,Y
make_sample <- function(null) sapply(1:50, function(i) Expm(mu[, , i], pdSpecEst:::T_coeff_inv(ifelse(null, 1, 0.5) * rexp(4) - 1, mu[, , i])), simplify = "array") 

X3 <- make_sample(null = T)
Y3 <- make_sample(null = T) ## null true
Y3a <- make_sample(null = F) ## null false (scale change)

## Signed-rank test
pdRankTests(array(c(X3, Y3), dim = c(2, 2, 100)), test = "signed.rank")[1:4] ## null true
pdRankTests(array(c(X3, Y3a), dim = c(2, 2, 100)), test = "signed.rank")[2] ## null false


## ------------------------------------------------------------------------
## Signed-rank test for equivalence of spectra
## ARMA(1,1) process: Example 11.4.1 in (Brockwell and Davis, 1991)
Phi <- array(c(0.7, 0, 0, 0.6, rep(0, 4)), dim = c(2, 2, 2))
Theta <- array(c(0.5, -0.7, 0.6, 0.8, rep(0, 4)), dim = c(2, 2, 2))
Sigma <- matrix(c(1, 0.71, 0.71, 2), nrow = 2)
pgram <- function(Sigma) pdPgram(rARMA(2^9, 2, Phi, Theta, Sigma)$X)$P ## HPD periodogram

## Null is true
pdRankTests(array(c(pgram(Sigma), pgram(Sigma)), dim = c(2, 2, 2^8)), test = "signed.rank")[2]

## Null is false
pdRankTests(array(c(pgram(Sigma), pgram(0.5 * Sigma)), dim = c(2, 2, 2^8)), test = "signed.rank")[2]

## ------------------------------------------------------------------------
## Null is true
data3 <- replicate(200, Expm(diag(2), pdSpecEst:::E_coeff_inv(rnorm(4)))) ## pointwise samples
data4 <- replicate(100, sapply(1:5, function(i) Expm(i * diag(2), pdSpecEst:::T_coeff_inv(rnorm(4), i * diag(2))), simplify = "array")) ## curve samples

## Null is false
data3a <- sapply(1:200, function(j) Expm(diag(2), pdSpecEst:::E_coeff_inv(((200 - j) / 200 + j * 2 / 200) * rnorm(4))), simplify = "array") ## pointwise trend in scale
data4a <- sapply(1:100, function(j) sapply(1:5, function(i) Expm(i * diag(2), pdSpecEst:::T_coeff_inv(((100 - j) / 100 + j * 2 / 100) * rnorm(4), i * diag(2))), simplify = "array"), simplify = "array") ## curve trend in scale

## Bartels-von Neumann test
pdRankTests(data3, test = "bartels")[1:4] ## null true (pointwise)
pdRankTests(data4, test = "bartels")[2] ## null true (curve)
pdRankTests(data3a, test = "bartels")[2] ## null false (pointwise)
pdRankTests(data4a, test = "bartels")[2] ## null false (curve)

## ---- echo = FALSE-------------------------------------------------------
d <- round(exp(seq(from = log(2), to = log(32), length = 5)), digits = 0)
n <- round(exp(seq(from = log(50), to = log(500), length = 5)), digits = 0)
par(mfrow = c(1,2), mgp = c(2,0.5,0), mar = c(4,3,2,2))
plot(range(d), range(pdSpecEst:::rs_time[1:3,,], na.rm=T), type="n", log = "xy", main = "Increasing dimension d x d, (n = 100)", xlab = "d", ylab = "Time (ms)")
for(i in 1:3){
  lines(d, pdSpecEst:::rs_time[i,,1], col = i, lty = i)
  points(d, pdSpecEst:::rs_time[i,,1], col = i, pch = 15+i)
}
legend("bottomright", legend = c("gdd", "zonoid", "spatial"), lty = 1:3, pch = 16:18, col = 1:3, cex = 1, bty = "n", y.intersp = 1)
plot(range(n), range(pdSpecEst:::rs_time[1:3,,], na.rm=T), type="n", log = "xy", main = "Increasing sample size n, (d = 4)", xlab = "n", ylab = "Time (ms)")
for(i in 1:3){
  lines(n, pdSpecEst:::rs_time[i,,2], col = i, lty = i)
  points(n, pdSpecEst:::rs_time[i,,2], col = i, pch = 15+i)
}
legend("bottomright", legend = c("gdd", "zonoid", "spatial"), lty = 1:3, pch = 16:18, col = 1:3, cex = 1, bty = "n", y.intersp = 1)

## ---- echo = FALSE-------------------------------------------------------
d <- round(exp(seq(from = log(2), to = log(32), length = 5)), digits = 0)
n <- round(exp(seq(from = log(25), to = log(250), length = 5)), digits = 0)
par(mfrow = c(1,2), mgp = c(2,0.5,0), mar = c(4,3,2,2))
plot(range(d), range(pdSpecEst:::kw_time[1:3,,], na.rm=T), type="n", log = "xy", main = "Increasing dimension d x d, (n = 50)", xlab = "d", ylab = "Time (ms)")
for(i in 1:3){
  lines(d, pdSpecEst:::kw_time[i,,1], col = i, lty = i)
  points(d, pdSpecEst:::kw_time[i,,1], col = i, pch = 15+i)
}
legend("bottomright", legend = c("gdd", "zonoid", "spatial"), lty = 1:3, pch = 16:18, col = 1:3, cex = 1, bty = "n", y.intersp = 1)
plot(range(n), range(pdSpecEst:::kw_time[1:3,,], na.rm=T), type="n", log = "xy", main = "Increasing sample size n, (d = 4)", xlab = "n", ylab = "Time (ms)")
for(i in 1:3){
  lines(n, pdSpecEst:::kw_time[i,,2], col = i, lty = i)
  points(n, pdSpecEst:::kw_time[i,,2], col = i, pch = 15+i)
}
legend("bottomright", legend = c("gdd", "zonoid", "spatial"), lty = 1:3, pch = 16:18, col = 1:3, cex = 1, bty = "n", y.intersp = 1)

## ---- echo = FALSE-------------------------------------------------------
d <- round(exp(seq(from = log(2), to = log(32), length = 5)), digits = 0)
n <- round(exp(seq(from = log(50), to = log(500), length = 5)), digits = 0)
par(mfrow = c(1,2), mgp = c(2,0.5,0), mar = c(4,3,2,2))
plot(range(d), range(pdSpecEst:::bart_time, na.rm=T), type="n", log = "xy", main = "Increasing dimension d x d, (n = 200)", xlab = "d", ylab = "Time (ms)")
for(i in 1:3){
  lines(d, pdSpecEst:::bart_time[i,,1], col = i, lty = i)
  points(d, pdSpecEst:::bart_time[i,,1], col = i, pch = 15+i)
}
legend("bottomright", legend = c("gdd", "zonoid", "spatial"), lty = 1:3, pch = 16:18, col = 1:3, cex = 1, bty = "n", y.intersp = 1)
plot(range(n), range(pdSpecEst:::bart_time, na.rm=T), type="n", log = "xy", main = "Increasing sample size n, (d = 4)", xlab = "n", ylab = "Time (ms)")
for(i in 1:3){
  lines(n, pdSpecEst:::bart_time[i,,2], col = i, lty = i)
  points(n, pdSpecEst:::bart_time[i,,2], col = i, pch = 15+i)
}
legend("bottomright", legend = c("gdd", "zonoid", "spatial"), lty = 1:3, pch = 16:18, col = 1:3, cex = 1, bty = "n", y.intersp = 1)

## ---- echo = FALSE-------------------------------------------------------
d <- round(exp(seq(from = log(2), to = log(128), length = 10)), digits = 0)
n <- round(exp(seq(from = log(50), to = log(500), length = 10)), digits = 0)
par(mfrow = c(1,2), mgp = c(2,0.5,0), mar = c(4,3,2,2))
plot(range(d), range(pdSpecEst:::sr_time, na.rm=T), type="n", log = "xy", main = "Increasing dimension d x d, (n = 100)", xlab = "d", ylab = "Time (ms)")
  lines(d, pdSpecEst:::sr_time[,1], col = 1, lty = 1)
  points(d, pdSpecEst:::sr_time[,1], col = 1, pch = 16)
plot(range(n), range(pdSpecEst:::sr_time, na.rm=T), type="n", log = "xy", main = "Increasing sample size n, (d = 6)", xlab = "n", ylab = "Time (ms)")
  lines(n, pdSpecEst:::sr_time[,2], col = 2, lty = 2)
  points(n, pdSpecEst:::sr_time[,2], col = 2, pch = 17)

