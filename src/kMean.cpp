# define ARMA_DONT_PRINT_ERRORS
# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_mat kMean(arma::cx_mat M, arma::vec mu) {

  int d = M.n_cols;

  int n = (int)M.n_rows / d;

  arma::cx_mat M1 = M.head_rows(d);

  arma::cx_mat Mi(d, d);

  arma::cx_mat M1sq(d, d);

  arma::cx_mat M1isq(d, d);

  arma::vec eigval;

  arma::cx_mat eigvec;

  double mu1;

  for(int i = 1; i < n; ++i) {

    if (i % 10000 == 0) {

      Rcpp::checkUserInterrupt();

    }

    Mi = M.rows(i * d, (i + 1) * d - 1);

    M1sq = arma::sqrtmat_sympd(M1);

    M1isq = arma::inv_sympd(M1sq);

    mu1 = arma::sum(mu.head(i + 1));

    arma::eig_sym(eigval, eigvec, M1isq * Mi * M1isq);

    arma::cx_mat M11 = M1sq * eigvec;

    if(mu1 == 0){

      M1 = M11 * arma::diagmat(arma::pow(eigval, mu[i])) * M11.t();

    } else {

      M1 = M11 * arma::diagmat(arma::pow(eigval, mu[i] / mu1)) * M11.t();

    }

  }

  return M1;
}
