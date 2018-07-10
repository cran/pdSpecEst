#define ARMA_DONT_PRINT_ERRORS
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' Geodesic midpoint between HPD matrices
//'
//' \code{Mid} calculates the geodesic midpoint between two HPD matrices under the
//' affine-invariant Riemannian metric as in \insertCite{B09}{pdSpecEst}[Chapter 6].
//'
//' @param A,B Hermitian positive definite matrices (of equal dimension).
//'
//' @examples
//'  ## Generate two random HPD matrices
//'  a <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
//'  A <- t(Conj(a)) %*% a
//'  b <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
//'  B <- t(Conj(b)) %*% b
//'  ## Compute midpoint
//'  Mid(A, B)
//'  ## Midpoint coincides with two-point intrinsic Karcher mean
//'  all.equal(pdMean(array(c(A, B), dim = c(3, 3, 2))), Mid(A, B))
//'
//' @references
//' \insertAllCited{}
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

arma::cx_mat Sqrt(arma::cx_mat M) {
  return arma::sqrtmat_sympd(M);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_mat Chol_C(arma::cx_mat M, bool bias, bool inverse) {
  try{
    // Declarations
    int d = M.n_rows;
    arma::cx_mat M1(d, d);
    // Forward Cholesky
    if(!inverse) {
      M1 = arma::trans(arma::chol(M));
    }
    // Backward Cholesky w/o bias
    if(inverse && !bias) {
      M1 = M * arma::trans(M);
    }
    // Backward Cholesky w bias
    if(inverse && bias) {
      arma::vec bias_vec(d);
      for(int k = 0; k < d; ++k) {
        bias_vec(k) = std::tgamma((double)(d - k + 0.5)) / (std::sqrt((double)d) * std::tgamma((double)(d - k)));
      }
      M1(0, 0) = (M(0, 0) * M(0, 0)) / (bias_vec(0) * bias_vec(0));
      for(int k = 0; k < (d - 1); ++k) {
        M1(k + 1, arma::span(0, k)) = M(k + 1, arma::span(0, k)) *
          arma::inv_sympd(M.submat(0, 0, k, k)) * M1.submat(0, 0, k, k);
        M1(arma::span(0, k), k + 1) = arma::trans(M1(k + 1, arma::span(0, k)));
        M1(k + 1, k + 1) = (M(k + 1, k + 1) * M(k + 1, k + 1)) / (bias_vec(k + 1) * bias_vec(k + 1)) +
          arma::as_scalar(M1(k + 1, arma::span(0, k)) * arma::inv_sympd(M1.submat(0, 0, k, k)) *
          arma::trans(M1(k + 1, arma::span(0, k))));
      }
    }
    return M1;

    // catch errors
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    Rcpp::exception("c++ exception (unknown reason)");
  }
  return arma::zeros<arma::cx_mat>(1,1);
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

arma::cx_cube wavPyr_C(arma::cx_cube P, int L, int J, arma::ivec Nj, std::string metric) {
  try{
    // Set parameters
    int d = P.n_rows;
    int n = P.n_slices;
    int N = (2 * L + 1) * n;
    arma::cx_cube P1 = P;
    arma::cx_cube P_rev(d, d, n);
    arma::cx_cube P_per(d, d, N);
    arma::cx_cube M_per(d, d, arma::sum(Nj));

    // Transform according to metric
    for(int i = 0; i < n; ++i) {
      if(metric == "logEuclidean") {
        P1.slice(i) = arma::logmat_sympd(P.slice(i));
      }
      else if(metric == "Cholesky") {
        P1.slice(i) = arma::trans(arma::chol(P.slice(i)));
      }
      else if(metric == "rootEuclidean") {
        P1.slice(i) = arma::sqrtmat_sympd(P.slice(i));
      }
    }
    // Construct reverse cube
    P_rev.slice(0) = P1.slice(n - 1);
    for(int i = 1; i < n; ++i) {
      P_rev.slice(i) = P1.slice(n - i);
    }
    // Construct periodized cube
    for(int l = 0; l < (2 * L + 1); ++l) {
      if((L % 2 == 0) && (l % 2 == 0)) {
        P_per.slices(l * n, (l + 1) * n - 1) = P1;
      }
      else if((L % 2 == 0) && (l % 2 != 0)) {
        P_per.slices(l * n, (l + 1) * n - 1) = P_rev;
      }
      else if((L % 2 != 0) && (l % 2 == 0)) {
        P_per.slices(l * n, (l + 1) * n - 1) = P_rev;
      }
      else if((L % 2 != 0) && (l % 2 != 0)) {
        P_per.slices(l * n, (l + 1) * n - 1) = P1;
      }
    }

    M_per.slices(0, N - 1) = P_per;
    for(int j = 1; j <= J; ++j) {
      int len = arma::sum(Nj.head(j));
      for(int k = 0; k < Nj(j); ++k) {
        if(metric == "Riemannian") {
          M_per.slice(len + k) = Mid(M_per.slice(len - Nj(j - 1) + 2 * k),
                      M_per.slice(len - Nj(j - 1) + 2 * k + 1));
        }
        else {
          M_per.slice(len + k) = (M_per.slice(len - Nj(j - 1) + 2 * k) +
            M_per.slice(len - Nj(j - 1) + 2 * k + 1)) / 2;
        }
      }
    }
    return M_per;

    // catch errors
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    Rcpp::exception("c++ exception (unknown reason)");
  }
  return arma::zeros<arma::cx_cube>(1,1,1);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_cube wavCoeff_C(arma::cx_cube M1, arma::cx_cube M0, double j, std::string metric) {
  try{
    // Set parameters
    int d = M1.n_cols;
    int n = M1.n_slices;
    arma::cx_cube W(d, d, 2 * n);
    arma::cx_mat Msqrt(d, d);
    arma::cx_mat Misqrt(d, d);

    for(int k = 0; k < n; ++k) {
      // Compute (non-)whitened wavelet coeff's
      if(metric == "Riemannian") {
        // Riemannian coefficients
        Msqrt = arma::sqrtmat_sympd(M1.slice(k));
        Misqrt = arma::inv_sympd(Msqrt);
        W.slice(k) = std::pow((double)2, (double)(-j / 2)) * arma::logmat_sympd(Misqrt * M0.slice(k) * Misqrt);
        W.slice(k + n) = Msqrt * W.slice(k) * Msqrt;
      }
      else {
        // Euclidean coefficients
        W.slice(k) = std::pow((double)2, (double)(-j / 2)) * (M0.slice(k) - M1.slice(k));
        W.slice(k + n) = W.slice(k);
      }
    }
    return W;

    // catch errors
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    Rcpp::exception("c++ exception (unknown reason)");
  }
  return arma::zeros<arma::cx_cube>(1,1,1);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_cube reconstr_C(arma::cx_cube M1, arma::cx_cube M0, arma::cx_cube Dj, double j,
                         int nj, bool in_sample, int L1, std::string metric) {
  try{
    // Initialize variables
    int d = M1.n_cols;
    arma::cx_cube M2(d, d, 2 * nj);
    arma::cx_mat Msqrt(d, d);
    arma::cx_mat Misqrt(d, d);

    for(int k = 0; k < nj; ++k) {
      // Reconstruct odd midpoints from non-zero wav. coeffs
      if((arma::norm(Dj.slice(k), "inf") > 1E-10) && in_sample) {
        if(metric == "Riemannian") {
          Msqrt = arma::sqrtmat_sympd(M1.slice(2 * k + 1));
          Misqrt = arma::inv_sympd(Msqrt);
          M2.slice(2 * k + 1) = Msqrt * arma::expmat_sym(std::pow((double)2, (double)(j/2)) *
            Misqrt * Dj.slice(k) * Misqrt) * Msqrt;
        }
        else {
          M2.slice(2 * k + 1) = std::pow((double)2, (double)(j/2)) * Dj.slice(k) + M1.slice(2 * k + 1);
        }
      }
      else {
        // Reconstruct odd midpoints from zero wav.coeffs
        M2.slice(2 * k + 1) = M1.slice(2 * k + 1);
      }
      // Reconstruct even midpoints from odd midpoints + coarse midpoints
      if(metric == "Riemannian") {
        M2.slice(2 * k) = M0.slice(k + L1) * arma::inv_sympd(M2.slice(2 * k + 1)) * M0.slice(k + L1);
      }
      else {
        M2.slice(2 * k) = 2 * M0.slice(k + L1) - M2.slice(2 * k + 1);
      }
    }
    return M2;

    // catch errors
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    Rcpp::exception("c++ exception (unknown reason)");
  }
  return arma::zeros<arma::cx_cube>(1,1,1);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_cube reconstr2D_C(arma::cx_cube M1, arma::cx_cube D, double j,
                           arma::ivec n, bool in_sample, std::string metric) {
  try{
    // Initialize variables
    int d = M1.n_cols;
    arma::cx_cube M2(d, d, n(0) * n(1));
    arma::cx_mat Msqrt(d, d);
    arma::cx_mat Misqrt(d, d);

    for(int k1 = 0; k1 < n(0); ++k1) {
      for(int k2 = 0; k2 < n(1); ++k2) {
        // Reconstruct midpoints from non-zero wav. coeffs
        if((arma::norm(D.slice(k2 * n(0) + k1), "inf") > 1E-10) && in_sample) {
          if(metric == "Riemannian") {
            // Riemannian metric
            Msqrt = arma::sqrtmat_sympd(M1.slice(k2 * n(0) + k1));
            Misqrt = arma::inv_sympd(Msqrt);
            M2.slice(k2 * n(0) + k1) = Msqrt * arma::expmat_sym(std::pow((double)2, (double)(j/2)) *
              Misqrt * D.slice(k2 * n(0) + k1) * Misqrt) * Msqrt;
          }
          else {
            // Euclidean metric
            M2.slice(k2 * n(0) + k1) = std::pow((double)2, (double)(j/2)) * D.slice(k2 * n(0) + k1) +
              M1.slice(k2 * n(0) + k1);
          }
        }
        else {
          // Reconstruct midpoints from zero wav. coeffs
          M2.slice(k2 * n(0) + k1) = M1.slice(k2 * n(0) + k1);
        }
      }
    }
    return M2;

    // catch errors
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    Rcpp::exception("c++ exception (unknown reason)");
  }
  return arma::zeros<arma::cx_cube>(1,1,1);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_cube Ptransf2D_C(arma::cx_cube P, bool inverse, bool chol_bias, std::string metric) {
  try{
    // Set parameters
    int d = P.n_rows;
    int n = P.n_slices;
    arma::cx_cube P1(d, d, n);

    if(!inverse) {
      // Transform according to metric
      for(int i = 0; i < n; ++i) {
        if(metric == "logEuclidean") {
          P1.slice(i) = arma::logmat_sympd(P.slice(i));
        }
        else if(metric == "Cholesky") {
          P1.slice(i) = arma::trans(arma::chol(P.slice(i)));
        }
        else if(metric == "rootEuclidean") {
          P1.slice(i) = arma::sqrtmat_sympd(P.slice(i));
        }
      }
    } else {
      // Transform back according to metric
      for(int i = 0; i < n; ++i) {
        if(metric == "logEuclidean") {
          P1.slice(i) = arma::expmat_sym(P.slice(i));
        }
        else if(metric == "Cholesky") {
          P1.slice(i) = Chol_C(P.slice(i), chol_bias, true);
        }
        else if(metric == "rootEuclidean") {
          P1.slice(i) = arma::trans(P.slice(i)) * P.slice(i);
        }
      }
    }
    return P1;

    // catch errors
  } catch(std::exception &ex) {
    forward_exception_to_r(ex);
  } catch(...) {
    Rcpp::exception("c++ exception (unknown reason)");
  }
  return arma::zeros<arma::cx_cube>(1,1,1);
}

