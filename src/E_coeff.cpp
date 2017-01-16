# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::vec E_coeff(arma::cx_mat H, arma::cx_cube E) {

 int d = H.n_rows;

 arma::vec coeff(d * d);

 for(int i=0; i < (d * d); i++){

   arma::cx_mat Ei = E.slice(i);

   arma::cx_double coeff0 = accu(H % conj(Ei));

   coeff[i] = real(coeff0);

 }

 return coeff;

}



