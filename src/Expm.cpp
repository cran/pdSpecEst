# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' Exponential map
//'
//' \code{Expm(P, H)} computes the projection of a Hermitian matrix \code{H} from the tangent space at a Hermitian
//' PD matrix \code{P} to the Riemannian manifold of Hermitian PD matrices via the
//' exponential map as in (Pennec, 2006). This is the unique inverse of the logarithmic map \code{\link{Logm}}.
//'
//' @param P a Hermitian positive definite matrix.
//' @param H a Hermitian matrix (of equal dimension as \code{P}).
//'
//' @examples
//'  H <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
//'  diag(H) <- rnorm(3)
//'  H[lower.tri(H)] <- t(Conj(H))[lower.tri(H)]
//'  p <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
//'  P <- t(Conj(p)) %*% p
//'  Expm(P, H)
//'
//' @references
//' Pennec, X. (2006). Intrinsic statistics on Riemannian manifolds: Basic tools for geometric
//' measurements. \emph{Journal of Mathematical Imaging and Vision} 25(1), 127-154.
//'
//' @seealso \code{\link{Logm}}
//'
//' @export
// [[Rcpp::export()]]

arma::cx_mat Expm(arma::cx_mat P, arma::cx_mat H) {

  arma::cx_mat P1 = arma::sqrtmat_sympd(P);

  arma::cx_mat P2 = arma::inv_sympd(P1);

  arma::cx_mat P3 = arma::expmat_sym(P2 * H * P2);

  return P1 * P3 * P1;

}


