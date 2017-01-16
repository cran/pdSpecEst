# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_mat Sqrt(arma::cx_mat M) {

  return arma::sqrtmat_sympd(M);

}
