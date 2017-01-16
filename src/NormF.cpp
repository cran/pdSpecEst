# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

double NormF(arma::cx_mat M) {

  return norm(M, "fro");

}


