# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

double RiemmDist(arma::cx_mat A, arma::cx_mat B) {

  arma::cx_mat A1 = arma::sqrtmat_sympd(A);

  arma::cx_mat A2 = arma::inv_sympd(A1);

  arma::cx_mat A3 = arma::logmat_sympd(A2 * B * A2);

  return norm(A3, "fro");

}
