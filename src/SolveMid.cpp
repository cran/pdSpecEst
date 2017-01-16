# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_mat solveMid(arma::cx_mat B, arma::cx_mat C) {

  return C * arma::inv_sympd(B) * C;

}
