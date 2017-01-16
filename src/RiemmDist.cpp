# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' Riemannian distance HPD-matrices
//'
//' \code{RiemmDist} calculates the natural Riemannian distance between two Hermitian PD matrices as
//' in (Bhatia, 2009, Chapter 6).
//'
//' @param A,B Hermitian positive-definite matrices (of equal dimension).
//'
//' @examples
//'  a <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
//'  A <- t(Conj(a)) %*% a
//'  b <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
//'  B <- t(Conj(b)) %*% b
//'  RiemmDist(A, B)
//'
//' @references Bhatia, R. (2009). \emph{Positive Definite Matrices}. New Jersey: Princeton University Press.
//'
//' @export
// [[Rcpp::export()]]

double RiemmDist(arma::cx_mat A, arma::cx_mat B) {

  arma::cx_mat A1 = arma::sqrtmat_sympd(A);

  arma::cx_mat A2 = arma::inv_sympd(A1);

  arma::cx_mat A3 = arma::logmat_sympd(A2 * B * A2);

  return norm(A3, "fro");

}
