# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::mat ARMA(arma::cube Phi, arma::cube Theta, arma::mat Z, int len) {

    int d = Phi.n_rows;

    arma::mat X(d, len);

    X.col(0) = arma::zeros<arma::vec>(d);
    X.col(1) = arma::zeros<arma::vec>(d);

    for(int i=2; i < len; i++){

      X.col(i) = Phi.slice(0) * X.col(i-1) + Phi.slice(1) * X.col(i-2) +
                    Z.col(i) + Theta.slice(0) * Z.col(i-1) + Theta.slice(1) * Z.col(i-2);

    }

    return X;
}

