# define ARMA_DONT_PRINT_ERRORS
# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::mat ARMA(arma::cube Phi, arma::cube Theta, arma::mat Z, int len) {
  // Simulate arma process
  try{
    int d = Phi.n_rows;
    arma::mat X(d, len);
    X.col(0) = arma::zeros<arma::vec>(d);
    X.col(1) = arma::zeros<arma::vec>(d);
    for(int i=2; i < len; ++i){
      X.col(i) = Phi.slice(0) * X.col(i-1) + Phi.slice(1) * X.col(i-2) +
        Z.col(i) + Theta.slice(0) * Z.col(i-1) + Theta.slice(1) * Z.col(i-2);
    }
    return X;
    // catch errors
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    Rcpp::exception("c++ exception (unknown reason)");
  }
  return arma::zeros<arma::mat>(1, 1); // not reached
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_cube pgram_C(arma::cx_mat X, int B, arma::cx_mat h, std::string method, bool is_2D) {
  // Compute fast periodogram
  try{
    int d = X.n_cols;
    int n = X.n_rows;
    if(method == "bartlett" && !is_2D){
      // Bartlett averaged periodogram
      int n_B = n / B;
      int m_B = n_B / 2;
      arma::cx_mat dft(n_B, d);
      arma::cx_mat dft_t(n_B, d);
      arma::cx_cube per = arma::zeros<arma::cx_cube>(d, d, m_B);
      arma::cx_cube per_b(d, d, m_B);
      for(int b = 0; b < B; ++b) {
        dft = arma::fft(X.rows(n_B * b, n_B * (b + 1) - 1)) / std::sqrt((double)n_B);
        dft_t = dft.t();
        for(int i = 0; i < m_B; ++i) {
          per_b.slice(i) = dft_t.col(i) * dft.row(i);
        }
        per += per_b / (2 * B * arma::datum::pi);
      }
      return arma::conj(per);
    }
    else if(method == "multitaper"){
      // DPSS multitaper periodogram
      int m;
      if(!is_2D) {
        m = n / 2;
      }
      else {
        m = n;
      }
      arma::cx_mat dft(m, d);
      arma::cx_mat dft_t(m, d);
      arma::cx_cube per = arma::zeros<arma::cx_cube>(d, d, m);
      arma::cx_cube per_b(d, d, m);
      for(int b = 0; b < B; ++b) {
        arma::cx_mat X_b = X;
        for(int k = 0; k < d; ++k) {
          X_b.col(k) %= h.col(b);
        }
        dft = arma::fft(X_b) / std::sqrt((double)n);
        dft_t = dft.t();
        for(int i = 0; i < m; ++i) {
          per_b.slice(i) = dft_t.col(i) * dft.row(i);
        }
        per += per_b / (2 * B * arma::datum::pi);
      }
      return arma::conj(per);
    }
    // catch errors
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    Rcpp::exception("c++ exception (unknown reason)");
  }
  return arma::zeros<arma::cx_cube>(1,1,1);
}
