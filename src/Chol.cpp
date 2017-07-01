# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_mat Chol(arma::cx_mat M) {

  return arma::chol(M);

}
