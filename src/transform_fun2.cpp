# define ARMA_DONT_PRINT_ERRORS
# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' Riemannian HPD exponential map
//'
//' \code{Expm(P, H)} computes the projection of a Hermitian matrix \code{H} from the tangent space at a Hermitian
//' PD matrix \code{P} to the manifold of Hermitian PD matrices equipped with the affine-invariant Riemannian metric
//' via the exponential map as in e.g., \insertCite{PFA05}{pdSpecEst}. This is the unique inverse of the Riemannian
//' logarithmic map \code{\link{Logm}}.
//'
//' @param P a Hermitian positive definite matrix.
//' @param H a Hermitian matrix (of equal dimension as \code{P}).
//'
//' @examples
//'  ## Generate random Hermitian matrix
//'  H <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
//'  diag(H) <- rnorm(3)
//'  H[lower.tri(H)] <- t(Conj(H))[lower.tri(H)]
//'  ## Generate random HPD matrix
//'  p <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
//'  P <- t(Conj(p)) %*% p
//'  ## Compute exponential map
//'  Expm(P, H)
//'
//' @references
//' \insertAllCited{}
//'
//' @seealso \code{\link{Logm}, \link{pdParTrans}}
//'
//' @export
// [[Rcpp::export()]]

arma::cx_mat Expm(arma::cx_mat P, arma::cx_mat H) {
  // Set parameters
  int d = P.n_cols;
  // Compute exponential map
  if(arma::norm(P - arma::eye<arma::mat>(d, d), "inf") < 1E-10) {
    return arma::expmat_sym(H);
  } else {
    arma::cx_mat P1 = arma::sqrtmat_sympd(P);
    arma::cx_mat P2 = arma::inv_sympd(P1);
    arma::cx_mat P3 = arma::expmat_sym(P2 * H * P2);
    return P1 * P3 * P1;
  }
}

// [[Rcpp::depends(RcppArmadillo)]]

//' Riemannian HPD logarithmic map
//'
//' \code{Logm(P, Q)} computes the projection of a Hermitian PD matrix \code{Q} in the manifold of HPD matrices
//' equipped with the affine-invariant Riemannian metric to the tangent space attached at the Hermitian PD matrix
//' \code{P} via the logarithmic map as in e.g., \insertCite{PFA05}{pdSpecEst}. This is the unique inverse of
//' the exponential map \code{\link{Expm}}.
//'
//' @param P a Hermitian positive definite matrix.
//' @param Q a Hermitian positive definite matrix (of equal dimension as \code{P}).
//'
//' @examples
//'  ## Generate two random HPD matrices
//'  q <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
//'  Q <- t(Conj(q)) %*% q
//'  p <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
//'  P <- t(Conj(p)) %*% p
//'  ## Compute logarithmic map
//'  Logm(P, Q)
//'
//' @references
//' Pennec, X. (2006). Intrinsic statistics on Riemannian manifolds: Basic tools for geometric
//' measurements. \emph{Journal of Mathematical Imaging and Vision} 25(1), 127-154.
//'
//' @seealso \code{\link{Expm}, \link{pdParTrans}}
//'
//' @export
// [[Rcpp::export()]]

arma::cx_mat Logm(arma::cx_mat P, arma::cx_mat Q) {
  // Set parameters
  int d = P.n_cols;
  // Compute logarithmic map
  if(arma::norm(P - arma::eye<arma::mat>(d, d), "inf") < 1E-10) {
    return arma::logmat_sympd(Q);
  } else {
    arma::cx_mat P1 = arma::sqrtmat_sympd(P);
    arma::cx_mat P2 = arma::inv_sympd(P1);
    arma::cx_mat P3 = arma::logmat_sympd(P2 * Q * P2);
    return P1 * P3 * P1;
  }
}

// [[Rcpp::depends(RcppArmadillo)]]

//' Riemannian HPD parallel transport
//'
//' \code{pdParTrans} computes the parallel transport on the manifold of HPD matrices
//' equipped with the affine-invariant Riemannian metric as described in e.g., Chapter 2 of \insertCite{C18}{pdSpecEst}. That is,
//' the function computes the parallel transport of a Hermitian matrix \code{W} in the tangent space
//' at the HPD matrix \code{P} along a geodesic curve in the direction of the Hermitian matrix \code{V}
//' in the tangent space at \code{P} for a unit time step.
//'
//' @param P a \eqn{(d,d)}-dimensional HPD matrix.
//' @param V a \eqn{(d,d)}-dimensional Hermitian matrix corresponding to a vector in the tangent space of \code{P}.
//' @param W a \eqn{(d,d)}-dimensional Hermitian matrix corresponding to a vector in the tangent space of \code{P}.
//'
//' @return a \eqn{(d,d)}-dimensional Hermitian matrix corresponding to the parallel transportation of \code{W} in
//' the direction of \code{V} along a geodesic curve for a unit time step.
//'
//' @examples
//' ## Transport the vector W to the tangent space at the identity
//' W <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
//' diag(W) <- rnorm(3)
//' W[lower.tri(W)] <- t(Conj(W))[lower.tri(W)]
//' p <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
//' P <- t(Conj(p)) %*% p
//'
//' pdParTrans(P, Logm(P, diag(3)), W) ## whitening transport
//'
//' @references
//' \insertAllCited{}
//'
//' @seealso \code{\link{Expm}, \link{Logm}}
//'
//' @export
// [[Rcpp::export()]]

arma::cx_mat pdParTrans(arma::cx_mat P, arma::cx_mat V, arma::cx_mat W) {
  // Compute parallel transport
  arma::cx_mat P1 = arma::sqrtmat_sympd(P);
  arma::cx_mat P2 = arma::inv_sympd(P1);
  arma::cx_mat P3 = P1 * arma::expmat_sym(0.5 * P2 * V * P2) * P2;
  return P3 * W * P3.t();

}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::vec E_coeff(arma::cx_mat H) {
  // Set parameters
  int d = H.n_rows;
  arma::vec coeff(d * d);
  // Basis transform
  for(int i = 0; i < d; ++i) {
    for(int j = 0; j < d; ++j) {
      if(i == j) {
        coeff(i * d + j) = real(H(i,i ));
      } else if(i > j) {
        coeff(i * d + j) = 2 / std::sqrt((double)2) * real(H(i, j));
      } else{
        coeff(i * d + j) = 2 / std::sqrt((double)2) * imag(H(i, j));
      }
    }
  }
  return coeff;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::vec T_coeff(arma::cx_mat H, arma::cx_mat y){

  // Set parameters
  int d = H.n_rows;
  arma::vec coeff(d * d);
  arma::cx_mat Ei = arma::zeros<arma::cx_mat>(d, d);
  arma::cx_mat y_sqrt = arma::sqrtmat_sympd(y);
  // Basis transform base y
  for(int i = 0; i < d; ++i) {
    for(int j = 0; j < d; ++j) {
      if(i == j) {
        Ei(i, i) = 1;
        Ei = y_sqrt * Ei * y_sqrt;
      } else if(i > j) {
        Ei(i, j) = Ei(j, i) = arma::cx_double(1 / std::sqrt((double)2), 0);
        Ei = y_sqrt * Ei * y_sqrt;
      } else{
        Ei(i, j) = arma::cx_double(0, 1 / std::sqrt((double)2));
        Ei(j, i) = arma::cx_double(0, -1 / std::sqrt((double)2));
        Ei = y_sqrt * Ei * y_sqrt;
      }
      coeff(i * d + j) = real(accu(H % conj(Ei)));
      Ei.zeros();
    }
  }
  return coeff;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::cx_mat E_coeff_inv(arma::vec coeff) {

  //Set parameters
  int d = (int)std::sqrt((double)coeff.size());
  arma::cx_mat M = arma::zeros<arma::cx_mat>(d, d);
  arma::cx_mat Ei = arma::zeros<arma::cx_mat>(d, d);
  //Inverse basis transform base Id
  for(int i = 0; i < d; ++i) {
    for(int j = 0; j < d; ++j) {
      if(i == j) {
        Ei(i, i) = 1;
      } else if(i > j) {
        Ei(i, j) = Ei(j, i) = arma::cx_double(1 / std::sqrt((double)2), 0);
      } else{
        Ei(i, j) = arma::cx_double(0, 1 / std::sqrt((double)2));
        Ei(j, i) = arma::cx_double(0, -1 / std::sqrt((double)2));
      }
      M += coeff(i * d + j) * Ei;
      Ei.zeros();
    }
  }
  return M;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
arma::cx_mat T_coeff_inv(arma::vec coeff, arma::cx_mat y) {
  //Set parameters
  int d = y.n_rows;
  arma::cx_mat M = arma::zeros<arma::cx_mat>(d, d);
  arma::cx_mat Ei = arma::zeros<arma::cx_mat>(d, d);
  arma::cx_mat y_isqrt = arma::inv_sympd(arma::sqrtmat_sympd(y));
  // Inverse basis transform base y
  for(int i = 0; i < d; ++i) {
    for(int j = 0; j < d; ++j) {
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
      M += coeff(i * d + j) * Ei;
      Ei.zeros();
    }
  }
  return M;
}


