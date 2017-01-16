# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_mat iSqrt(arma::cx_mat M) {

  arma::cx_mat M1 = arma::sqrtmat_sympd(M);

  return arma::inv_sympd(M1);

}
