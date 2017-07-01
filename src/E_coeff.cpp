# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::vec E_coeff(arma::cx_mat H) {

  int d = H.n_rows;

  arma::vec coeff(d * d);

  for(int i = 0; i < d; i++) {
    for(int j = 0; j < d; j++) {

        if(i == j) {
          coeff[i * d + j] = real(H(i,i ));
        } else if(i > j) {
          coeff[i * d + j] = 2 / sqrt(2) * real(H(i, j));
        } else{
          coeff[i * d + j] = 2 / sqrt(2) * imag(H(i, j));
        }
    }
  }

  return coeff;

}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::vec T_coeff(arma::cx_mat H, arma::cx_mat y){

  int d = H.n_rows;

  arma::vec coeff(d * d);
  arma::cx_mat Ei = arma::zeros<arma::cx_mat>(d, d);
  arma::cx_mat y_sqrt = arma::sqrtmat_sympd(y);

  for(int i = 0; i < d; i++) {
    for(int j = 0; j < d; j++) {

      if(i == j) {

        Ei(i, i) = 1;
        Ei = y_sqrt * Ei * y_sqrt;

      } else if(i > j) {

        Ei(i, j) = Ei(j, i) = arma::cx_double(1 / sqrt(2), 0);
        Ei = y_sqrt * Ei * y_sqrt;

      } else{

        Ei(i, j) = arma::cx_double(0, 1 / sqrt(2));
        Ei(j, i) = arma::cx_double(0, -1 / sqrt(2));
        Ei = y_sqrt * Ei * y_sqrt;

      }

      coeff[i * d + j] = real(accu(H % conj(Ei)));

      Ei.zeros();
    }
  }

  return coeff;

}

