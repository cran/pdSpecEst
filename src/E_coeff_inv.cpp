# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::cx_mat E_coeff_inv(arma::vec coeff) {

  int d = (int)std::sqrt((double)coeff.size());

  arma::cx_mat M = arma::zeros<arma::cx_mat>(d, d);
  arma::cx_mat Ei = arma::zeros<arma::cx_mat>(d, d);

  for(int i = 0; i < d; i++) {
    for(int j = 0; j < d; j++) {

      if(i == j) {

        Ei(i, i) = 1;

      } else if(i > j) {

        Ei(i, j) = Ei(j, i) = arma::cx_double(1 / std::sqrt((double)2), 0);

      } else{

        Ei(i, j) = arma::cx_double(0, 1 / std::sqrt((double)2));

        Ei(j, i) = arma::cx_double(0, -1 / std::sqrt((double)2));

      }

      M += coeff[i * d + j] * Ei;

      Ei.zeros();

    }
  }

  return M;

}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::cx_mat T_coeff_inv(arma::vec coeff, arma::cx_mat y) {

  int d = y.n_rows;

  arma::cx_mat M = arma::zeros<arma::cx_mat>(d, d);
  arma::cx_mat Ei = arma::zeros<arma::cx_mat>(d, d);
  arma::cx_mat y_isqrt = arma::inv_sympd(arma::sqrtmat_sympd(y));

  for(int i = 0; i < d; i++) {
    for(int j = 0; j < d; j++) {

      if(i == j) {

        Ei(i, i) = 1;
        Ei = y_isqrt * Ei * y_isqrt;

      } else if(i > j) {

        Ei(i, j) = Ei(j, i) = arma::cx_double(1 / std::sqrt((double)2), 0);
        Ei = y_isqrt * Ei * y_isqrt;

      } else{

        Ei(i, j) = arma::cx_double(0, 1 / std::sqrt((double)2));
        Ei(j, i) = arma::cx_double(0, -1 / std::sqrt((double)2));
        Ei = y_isqrt * Ei * y_isqrt;

      }

      M += coeff[i * d + j] * Ei;

      Ei.zeros();

    }
  }

  return M;

}

