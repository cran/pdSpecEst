# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_mat E_coeff_inv(arma::vec coeff, arma::cx_cube E) {

 int d2 = coeff.size();

 int d = E.n_rows;

 arma::cx_mat M = arma::zeros<arma::cx_mat>(d, d);

 for(int i=0; i < d2; i++){

   arma::cx_mat Ei = E.slice(i);

   M += coeff[i] * Ei;

 }

 return M;

}


