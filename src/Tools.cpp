#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' Geodesic midpoint between HPD-matrices
//'
//' \code{Mid} calculates the geodesic midpoint between two Hermitian PD matrices as in
//' (Bhatia, 2009, Chapter 6).
//'
//' @param A,B Hermitian positive definite matrices (of equal dimension).
//'
//' @examples
//'  a <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
//'  A <- t(Conj(a)) %*% a
//'  b <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
//'  B <- t(Conj(b)) %*% b
//'  Mid(A, B)
//' @references Bhatia, R. (2009). \emph{Positive Definite Matrices}. New Jersey: Princeton University Press.
//'
//' @seealso \code{\link{pdMean}}
//'
//' @export
// [[Rcpp::export()]]

arma::cx_mat Mid(arma::cx_mat A, arma::cx_mat B) {

  arma::cx_mat A1 = arma::sqrtmat_sympd(A);

  arma::cx_mat A2 = arma::inv_sympd(A1);

  arma::cx_mat C = A2 * B * A2;

  return A1 * arma::sqrtmat_sympd(C) * A1;

}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_mat solveMid(arma::cx_mat B, arma::cx_mat C) {

  return C * arma::inv_sympd(B) * C;

}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_mat Sqrt(arma::cx_mat M) {

  return arma::sqrtmat_sympd(M);

}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_mat iSqrt(arma::cx_mat M) {

  arma::cx_mat M1 = arma::sqrtmat_sympd(M);

  return arma::inv_sympd(M1);

}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

double NormF(arma::cx_mat M) {

  return arma::norm(M, "fro");

}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

double RiemmDist(arma::cx_mat A, arma::cx_mat B) {

  arma::cx_mat A1 = arma::sqrtmat_sympd(A);

  arma::cx_mat A2 = arma::inv_sympd(A1);

  arma::cx_mat A3 = arma::logmat_sympd(A2 * B * A2);

  return arma::norm(A3, "fro");

}

