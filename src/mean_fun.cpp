# define ARMA_DONT_PRINT_ERRORS
# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_mat Mid_w(arma::cx_mat A, arma::cx_mat B, double w, std::string metric) {
	arma::cx_mat res(A.n_rows, A.n_cols);
	if(metric == "Riemannian") {
		arma::cx_mat B1 = arma::sqrtmat_sympd(B);
		arma::cx_mat B2 = arma::inv_sympd(B1);
		arma::vec eigval;
		arma::cx_mat eigvec;
		arma::eig_sym(eigval, eigvec, B2 * A * B2);
		res = B1 * (eigvec * arma::diagmat(arma::pow(eigval, w)) * eigvec.t()) * B1;
	}
	else {
		res = (1-w) * B + w * A;
	}
	return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_mat pdMean_C_approx(arma::cx_cube M, arma::vec mu) {
	try{
		// initialize params
		int d = M.n_cols;
		int n = mu.size();
		arma::cx_mat M1 = M.slice(0);
		arma::cx_mat Mi(d, d);
		arma::cx_mat M1sq(d, d);
		arma::cx_mat M1isq(d, d);
		arma::vec eigval;
		arma::cx_mat eigvec;
		double mu1;
		// compute approx mean
		for(int i = 1; i < n; ++i) {
			if (i % 100 == 0) {
				Rcpp::checkUserInterrupt();
			}
			Mi = M.slice(i);
			M1sq = arma::sqrtmat_sympd(M1);
			M1isq = arma::inv_sympd(M1sq);
			mu1 = arma::sum(mu.head(i + 1));
			arma::eig_sym(eigval, eigvec, M1isq * Mi * M1isq);
			arma::cx_mat M11 = M1sq * eigvec;
			if(mu1 == 0){
				M1 = M11 * arma::diagmat(arma::pow(eigval, mu(i))) * M11.t();
			} else {
				M1 = M11 * arma::diagmat(arma::pow(eigval, mu(i) / mu1)) * M11.t();
			}
		}
		return M1;

		// catch exceptions
	} catch(std::exception &ex) {
		forward_exception_to_r(ex);
	} catch(...) {
		Rcpp::exception("c++ exception (unknown reason)");
	}
	return arma::zeros<arma::cx_mat>(1,1); //not reached
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_mat pdMean_C(arma::cx_mat M0, arma::cx_cube M, arma::vec mu,
		int maxit, double reltol) {
	try{
		// initialize params
		int d = M0.n_cols;
		int n = mu.size();
		arma::cx_mat M1 = M0;
		arma::cx_mat M0sqrt;
		arma::cx_mat M0isqrt;
		arma::cx_mat M1log;
		double reltol_m = reltol + 1;
		int i = 0;
		// run gradient descent
		while(reltol_m > reltol && i < maxit) {
			if (i % 100 == 0) {
				Rcpp::checkUserInterrupt();
			}
			M0 = M1;
			M0sqrt = arma::sqrtmat_sympd(M0);
			M0isqrt = arma::inv_sympd(M0sqrt);
			M1log = arma::zeros<arma::cx_mat>(d, d);
			for(int j = 0; j < n; ++j) {
				M1log += mu(j) * arma::logmat_sympd(M0isqrt * M.slice(j) * M0isqrt);
				if(M1log.has_nan()) {
					Rcpp::stop("c++ matrix functions failed, matrix possibly not positive definite");
				}
			}
			M1 = M0sqrt * arma::expmat_sym(M1log) * M0sqrt;
			reltol_m = arma::norm(arma::logmat_sympd(M0isqrt * M1 * M0isqrt), "fro");
			++i;
		}
		return M1;

		// catch errors
	} catch(std::exception &ex) {
		forward_exception_to_r(ex);
	} catch(...) {
		Rcpp::exception("c++ exception (unknown reason)");
	}
	return arma::zeros<arma::cx_mat>(1,1); //not reached
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_mat pdMedian_C(arma::cx_mat M0, arma::cx_cube M, arma::vec mu,
		int maxit, double reltol) {
	try{
		// initialize params
		int d = M0.n_cols;
		int n = mu.size();
		arma::cx_mat M1 = M0;
		arma::cx_mat M0sqrt;
		arma::cx_mat M0isqrt;
		arma::cx_mat M1log;
		arma::cx_mat Mlog;
		double dist;
		double W;
		double reltol_m = reltol + 1;
		int i = 0;
		// run Weiszfeld algorithm (gradient descent)
		while(reltol_m > reltol && i < maxit) {
			if (i % 100 == 0) {
				Rcpp::checkUserInterrupt();
			}
			M0 = M1;
			M0sqrt = arma::sqrtmat_sympd(M0);
			M0isqrt = arma::inv_sympd(M0sqrt);
			M1log = arma::zeros<arma::cx_mat>(d, d);
			W = 0;
			for(int j = 0; j < n; ++j) {
				Mlog = arma::logmat_sympd(M0isqrt * M.slice(j) * M0isqrt);
				if(Mlog.has_nan()) {
					Rcpp::stop("c++ matrix functions failed, matrix possibly not positive definite");
				}
				dist = arma::norm(Mlog, "fro");
				if(dist > 1E-10) {
					M1log += mu(j) / dist * Mlog;
					W += mu(j) / dist;
				}
			}
			M1 = M0sqrt * arma::expmat_sym(M1log / W) * M0sqrt;
			reltol_m = arma::norm(arma::logmat_sympd(M0isqrt * M1 * M0isqrt), "fro");
			++i;
		}
		return M1;

		// catch errors
	} catch(std::exception &ex) {
		forward_exception_to_r(ex);
	} catch(...) {
		Rcpp::exception("c++ exception (unknown reason)");
	}
	return arma::zeros<arma::cx_mat>(1,1); //not reached
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_mat Euclid_Median_C(arma::cx_mat M0, arma::cx_cube M, arma::vec mu,
		int maxit, double reltol) {
	try{
		// initialize parameters
		int n = mu.size();
		arma::cx_mat M1 = M0;
		double dist;
		double W;
		double reltol_m = reltol + 1;
		int i = 0;
		// run Weiszfeld algorithm (gradient descent)
		while(reltol_m > reltol && i < maxit) {
			if (i % 100 == 0) {
				Rcpp::checkUserInterrupt();
			}
			M0 = M1;
			M1.zeros();
			W = 0;
			for(int j = 0; j < n; ++j) {
				dist = arma::norm(M.slice(j) - M0, "fro");
				if(dist > 1E-10) {
					M1 += mu(j) / dist * M.slice(j);
					W += mu(j) / dist;
				}
			}
			M1 = M1 / W;
			reltol_m = arma::norm(M1 - M0, "fro");
			++i;
		}
		return M1;

		// catch errors
	} catch(std::exception &ex) {
		forward_exception_to_r(ex);
	} catch(...) {
		Rcpp::exception("c++ exception (unknown reason)");
	}
	return arma::zeros<arma::cx_mat>(1,1); //not reached
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_cube pdNeville_C(arma::cx_cube P, arma::vec X, arma::vec x, std::string metric) {

	int n_P = X.size() - 1;
	int n_x = x.size();
	int d = P.n_cols;
	arma::cx_cube res(d, d, n_x);
	for(int j = 0; j < n_x; ++j) {
		arma::cx_cube p = P;
		for(int k = 0; k < n_P; ++k) {
			for(int i = 0; i < (n_P - k); ++i) {
				p.slice(i) = Mid_w(p.slice(i + 1), p.slice(i), (x(j) - X(i)) / (X(i + k + 1) - X(i)), metric);
			}
		}
		res.slice(j) = p.slice(0);
	}
	return res;
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cx_cube impute_C(arma::cx_cube M0, arma::mat W, int L, bool inverse,
		std::string metric, std::string method) {
	try{
		// Set parameters
		int n = M0.n_slices;
		int d = M0.n_rows;
		int N = 2 * L + 1;
		arma::vec w = arma::ones<arma::vec>(N);
		arma::cx_cube M(d, d, N);
		arma::cx_cube M1 = arma::zeros<arma::cx_cube>(d, d, 2 * n);
		arma::cx_mat Msqrt(d, d);
		arma::cx_mat Misqrt(d, d);

		if(L == 0 || n == 1) {
			// Trivial refinement
			for(int k = 0; k < n; ++k) {
				M1.slice(2 * k + 1) = M0.slice(k);
				if(inverse) {
					M1.slice(2 * k) = M0.slice(k);
				}
			}
		}
		else if (L <= 4){
			// Nontrivial refinement with filter weights
			for(int k = 0; k < n; ++k) {
				// Coarse-scale midpoints and filter weights
				if((k - L) < 0) {
					M = M0.slices(0, N - 1);
					w = arma::conv_to<arma::vec>::from(W.row(2 * k + 1));
				}
				else if((k + L) > (n - 1)) {
					M = M0.slices(n - N, n - 1);
					w = arma::conv_to<arma::vec>::from(W.row(2 * (N - (n - k)) + 1));
				}
				else {
					M = M0.slices(k - L, k + L);
					w = arma::conv_to<arma::vec>::from(W.row(N));
				}
				// Compute predicted midpoints
				if(metric == "Riemannian") {
					// Riemannian weighted mean (odd locations)
					M1.slice(2 * k + 1) = pdMean_C_approx(M, w);
					if(method != "fast") {
						M1.slice(2 * k + 1) = pdMean_C(M1.slice(2 * k + 1), M, w, 1000, 1E-10);
					}
					// Prediction even locations
					if(inverse) {
						M1.slice(2 * k) = M0.slice(k) * arma::inv_sympd(M1.slice(2 * k + 1)) * M0.slice(k);
					}
				} else if(metric == "Riemannian-Rahman"){
					// Riemannian-Rahman weighted mean (odd locations)
					arma::cx_mat M_w = arma::zeros<arma::cx_mat>(d, d);
					Msqrt = arma::sqrtmat_sympd(M0.slice(k));
					Misqrt = arma::inv_sympd(Msqrt);

					for(int i = 0; i < N; ++i){
                    	// Logarithmic weighted average in tangent space
                    	M_w += w(i) * Msqrt * arma::logmat_sympd(Misqrt * M.slice(i) * Misqrt) * Msqrt;
					}
					// Exponential map of M_w at M0.slice(k)
					M1.slice(2 * k + 1) = Msqrt * arma::expmat_sym(Misqrt * M_w * Misqrt) * Msqrt;

					// Prediction even locations
					if(inverse) {
						M1.slice(2 * k) = M0.slice(k) * arma::inv_sympd(M1.slice(2 * k + 1)) * M0.slice(k);
					}
				} else {
					// Euclidean weighted mean (odd locations)
					arma::cx_mat M_w = arma::zeros<arma::cx_mat>(d, d);
					for(int i = 0; i < N; ++i){
						M_w += w(i) * M.slice(i);
					}
					M1.slice(2 * k + 1) = M_w;
					// Prediction even locations
					if(inverse) {
						M1.slice(2 * k) = 2 * M0.slice(k) - M1.slice(2 * k + 1);
					}
				}
			}
		}
		else if(L > 4) {
			// Nontrivial refinement with Neville's algorithm
			arma::cx_cube M_bar(d, d, N);
			arma::cx_cube M_neville(d, d, 1);
			arma::vec N_seq = arma::linspace(1, N, N);

			for(int k = 0; k < n; ++k) {
				// Coarse-scale midpoints
				if((k - L) < 0) {
					M = M0.slices(0, N - 1);
				}
				else if((k + L) > (n - 1)) {
					M = M0.slices(n - N, n - 1);
				}
				else {
					M = M0.slices(k - L, k + L);
				}
				// Compute predicted midpoints (Neville's algorithm)
				arma::cx_mat M_w = arma::zeros<arma::cx_mat>(d, d);
				for(int i = 0; i < N; ++i) {
					// Cumulative intrinsic means
					if(metric == "Riemannian") {
						M_bar.slice(i) = pdMean_C_approx(M.head_slices(i + 1), w.head(i + 1));
						// if(method != "fast") {
						//   M_bar.slice(i) = pdMean_C(M_bar.slice(i), M.head_slices(i + 1), w.head(i + 1), 1000, 1E-10);
						//   }
					}
					else {
						M_bar.slice(i) = arma::mean(M.head_slices(i + 1), 2);
					}
				}
				if((k - L) < 0) {
					// Prediction at left boundary, uneven locations
					M = M0.slices(0, N - 1);
					arma::vec k_seq(1);
					k_seq(0) = k + 0.5;
					M_neville = pdNeville_C(M_bar, N_seq, k_seq, metric);
					M1.slice(2 * k + 1) = Mid_w(M_neville.slice(0), M_bar.slice(k), -(2 * k + 1), metric);
				}
				else if((k + L) > (n - 1)) {
					// Prediction at right boundary, uneven locations
					M = M0.slices(n - N, n - 1);
					int k1 = N - (n - k);
					arma::vec k1_seq(1);
					k1_seq(0) = k1 + 0.5;
					M_neville = pdNeville_C(M_bar, N_seq, k1_seq, metric);
					M1.slice(2 * k + 1) = Mid_w(M_neville.slice(0), M_bar.slice(k1), -(2 * k1 + 1), metric);
				}
				else {
					// Prediction away from boundary, uneven locations
					M = M0.slices(k - L, k + L);
					arma::vec D_seq(1);
					D_seq(0) = L + 0.5;
					M_neville = pdNeville_C(M_bar, N_seq, D_seq, metric);
					M1.slice(2 * k + 1) = Mid_w(M_neville.slice(0), M_bar.slice(L), -N, metric);
				}
				// Prediction at even locations
				if(inverse) {
					if(metric == "Riemannian") {
						M1.slice(2 * k) = M0.slice(k) * arma::inv_sympd(M1.slice(2 * k + 1)) * M0.slice(k);
					}
					else {
						M1.slice(2 * k) = 2 * M0.slice(k) - M1.slice(2 * k + 1);
					}
				}
			}
		}
		return M1;

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

arma::cx_cube wavPyr2D_C(arma::cx_cube P, int n1, int n2, std::string metric) {
	try{
		// Initialize variables
		int d = P.n_rows;
		int n = P.n_slices;
		int len = n / 4;
		int len1, len2;
		if(n1 == (int)1 || n2 == (int)1) {
			len = n / 2;
		}
		arma::cx_cube M(d, d, len);
		arma::cx_cube M2d(d, d, 4);
		arma::vec w = arma::ones<arma::vec>(4) / 4;

		// Local averages curve of HPD matrices
		if(n1 == (int)1 || n2 == (int)1) {
			for(int k = 0; k < len; ++k) {
				M.slice(k) = Mid_w(P.slice(2 * k), P.slice(2 * k + 1), 0.5, metric);
			}
		}
		else {
			// Local averages surface of HPD matrices
			len1 = n1 / 2;
			len2 = n2 / 2;
			for(int k1 = 0; k1 < len1; ++k1) {
				for(int k2 = 0; k2 < len2; ++k2) {

					M2d.slice(0) = P.slice(2 * k2 * n1 + 2 * k1);
					M2d.slice(1) = P.slice(2 * k2 * n1 + 2 * k1 + 1);
					M2d.slice(2) = P.slice((2 * k2 + 1) * n1 + 2 * k1);
					M2d.slice(3) = P.slice((2 * k2 + 1) * n1 + 2 * k1 + 1);

					if(metric == "Riemannian") {
						M.slice(k2 * len1 + k1) = pdMean_C_approx(M2d, w);
					}
					else {
						M.slice(k2 * len1 + k1) = arma::mean(M2d, 2);
					}
				}
			}
		}
		return M;

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

arma::cx_cube impute2D_C(arma::cx_cube M0, arma::field<arma::mat> W, int n1, int n2,
		arma::ivec L, std::string metric, std::string method) {
	try{
		// Set parameters
		int d = M0.n_rows;
		arma::field<arma::cx_mat> M1_field(2 * n1, 2 * n2);
		arma::ivec N = 2 * L + 1;
		arma::ivec nbrs_x(2);
		arma::ivec nbrs_y(2);
		arma::ivec L_k(2);
		arma::field<arma::cx_mat> M0_field(n1, n2);
		for(int i1 = 0; i1 < n1; ++i1) {
			for(int i2 = 0; i2 < n2; ++i2) {
				M0_field(i1, i2) = M0.slice(i2 * n1 + i1);
			}
		}
		// Prediction across all locations (k1, k2)
		for(int k1 = 0; k1 < n1; ++k1) {
			for(int k2 = 0; k2 < n2; ++k2) {

				// Default refinement order away from boundary
				L_k = L;
				// Reduce refinement order at boundary locations
				if(k1 < L(0)) {
					L_k(0) = k1;
				}
				if(k1 > (n1 - 1 - L(0))) {
					L_k(0) = std::abs((int)(k1 - (n1 - 1)));
				}
				if(k2 < L(1)) {
					L_k(1) = k2;
				}
				if(k2 > (n2 - 1 - L(1))) {
					L_k(1) = std::abs((int)(k2 - (n2 - 1)));
				}

				// Neighboring coarse midpoints
				nbrs_x(0) = k1 - L_k(0);
				nbrs_x(1) = k1 + L_k(0);
				nbrs_y(0) = k2 - L_k(1);
				nbrs_y(1) = k2 + L_k(1);

				if(L_k(0) < 1 && L_k(1) < 1) {
					M1_field(2 * k1, 2 * k2) = M0_field(k1, k2);
					M1_field(2 * k1 + 1, 2 * k2) = M0_field(k1, k2);
					M1_field(2 * k1, 2 * k2 + 1) = M0_field(k1, k2);
					M1_field(2 * k1 + 1, 2 * k2 + 1) = M0_field(k1, k2);
				}
				else {
					// Midpoint prediction with available weights
					// Select subfield of neighboring coarse midpoints
					arma::field<arma::cx_mat> M_field = M0_field.subfield(nbrs_x(0), nbrs_y(0), nbrs_x(1), nbrs_y(1));
					// Transform subfield to cube
					arma::cx_cube M_cube(d, d, M_field.n_rows * M_field.n_cols);
					for(arma::uword i1 = 0; i1 < M_field.n_rows; ++i1){
						for(arma::uword i2 = 0; i2 < M_field.n_cols; ++i2) {
							M_cube.slice(i2 * M_field.n_rows + i1) = M_field(i1, i2);
						}
					}
					// Refinement weights
					arma::mat w = W(L_k(1) * 5 + L_k(0));
					if(metric == "Riemannian") {
						// Approximate Riemannian weighted averages
						M1_field(2 * k1, 2 * k2) = pdMean_C_approx(M_cube, w.col(0));
						M1_field(2 * k1 + 1, 2 * k2) = pdMean_C(M0_field(k1, k2), M_cube, w.col(1), 15, 0.01);
						M1_field(2 * k1, 2 * k2 + 1) = pdMean_C(M0_field(k1, k2), M_cube, w.col(2), 15, 0.01);
						M1_field(2 * k1 + 1, 2 * k2 + 1) = pdMean_C_approx(M_cube, w.col(3));
						if(method != "fast") {
							// Accurate Riemannian weighted averages
							M1_field(2 * k1, 2 * k2) = pdMean_C(M1_field(2 * k1, 2 * k2),
									M_cube, w.col(0), 1000, 1E-10);
							M1_field(2 * k1 + 1, 2 * k2) = pdMean_C(M1_field(2 * k1 + 1, 2 * k2),
									M_cube, w.col(1), 1000, 1E-10);
							M1_field(2 * k1, 2 * k2 + 1) = pdMean_C(M1_field(2 * k1, 2 * k2 + 1),
									M_cube, w.col(2), 1000, 1E-10);
							M1_field(2 * k1 + 1, 2 * k2 + 1) = pdMean_C(M1_field(2 * k1 + 1, 2 * k2 + 1),
									M_cube, w.col(3), 1000, 1E-10);
						}
					}
					else {
						// Euclidean weighted averages
						arma::cx_mat M_w0 = arma::zeros<arma::cx_mat>(d, d);
						arma::cx_mat M_w1 = arma::zeros<arma::cx_mat>(d, d);
						arma::cx_mat M_w2 = arma::zeros<arma::cx_mat>(d, d);
						arma::cx_mat M_w3 = arma::zeros<arma::cx_mat>(d, d);
						for(arma::uword i = 0; i < M_cube.n_slices; ++i){
							M_w0 += w(i, 0) * M_cube.slice(i);
							M_w1 += w(i, 1) * M_cube.slice(i);
							M_w2 += w(i, 2) * M_cube.slice(i);
							M_w3 += w(i, 3) * M_cube.slice(i);
						}
						M1_field(2 * k1, 2 * k2) = M_w0;
						M1_field(2 * k1 + 1, 2 * k2) = M_w1;
						M1_field(2 * k1, 2 * k2 + 1) = M_w2;
						M1_field(2 * k1 + 1, 2 * k2 + 1) = M_w3;
					}
				}
			}
		}
		// Transform field of predicted midpoints to cube
		arma::cx_cube M1(d, d, 4 * n1 * n2);
		for(arma::uword k = 0; k < M1.n_slices; ++k) {
			M1.slice(k) = M1_field(k % (2 * n1), k / (2 * n1));
		}
		return M1;

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

double pdDist_C(arma::cx_mat A, arma::cx_mat B, std::string metric) {
	try{
		arma::cx_mat A2;
		// compute distance
		if(metric == "Riemannian") {
			arma::cx_mat A1 = arma::inv_sympd(arma::sqrtmat_sympd(A));
			A2 = arma::logmat_sympd(A1 * B * A1);
		}
		else if (metric == "logEuclidean") {
			A2 = arma::logmat_sympd(A) - arma::logmat_sympd(B);
		}
		else if (metric == "Cholesky") {
			A2 = arma::chol(A) - arma::chol(B);
		}
		else if (metric == "rootEuclidean") {
			A2 = arma::sqrtmat_sympd(A) - arma::sqrtmat_sympd(B);
		}
		else {
			A2 = A - B;
		}
		if(A2.has_nan()) {
			Rcpp::stop("c++ function logmat_sympd() failed, matrix possibly not positive definite");
		}
		return arma::norm(A2, "fro");

		// catch exceptions
	} catch(std::exception &ex) {
		forward_exception_to_r(ex);
	} catch(...) {
		Rcpp::exception("c++ exception (unknown reason)");
	}
	return NA_REAL; //not reached
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::vec gdd_C(arma::cx_cube X, arma::cx_cube y, std::string metric, int S, std::string method) {
	try{
		// Initialize parameters
		int n = X.n_slices / S;
		int res_size;
		if(method == "depth.y") {
			res_size = 1;
		}
		else {
			res_size = S;
		}
		arma::vec res(res_size);

		if(method == "depth.y") {
			// Calculate individual (integrated) depth y wrt X
			double dd = 0.0;
			for(int s = 0; s < S; ++s) {
				// Integrated depth
				for(int i = 0; i < n; ++i) {
					dd += pdDist_C(y.slice(i), X.slice(s * n + i), metric);
				}
			}
			res(0) = std::exp(-dd / (n * S));
		}
		else if(method == "depth.X") {
			// Calculate all (integrated) depths X wrt X
			arma::mat dist_mat = arma::zeros<arma::mat>(S, S);
			for(int s1 = 0; s1 < S; ++s1) {
				for(int s2 = 0; s2 < s1; ++s2) {
					// Integrated depth
					for(int i = 0; i < n; ++i) {
						dist_mat(s1, s2) += pdDist_C(X.slice(s1 * n + i), X.slice(s2 * n + i), metric) / n;
					}
					dist_mat(s2, s1) = dist_mat(s1, s2);
				}
			}
			res = arma::conv_to<arma::vec>::from(arma::exp(-arma::mean(dist_mat, 0)));
		}
		return res;

		// catch errors
	} catch(std::exception &ex) {
		forward_exception_to_r(ex);
	} catch(...) {
		Rcpp::exception("c++ exception (unknown reason)");
	}
	return arma::zeros<arma::vec>(1);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]

arma::cube cMeans_C(arma::cx_cube M, arma::cx_cube centroids, int S, int K, double m,
		double eps, int max_iter, std::string metric, arma::mat dist_weights) {
	try{
		// Initialize parameters
		int d = M.n_rows;
		int n = M.n_slices / S;
		arma::mat dist_mat(S, K);
		arma::mat clust_mat(S, K);
		arma::cx_cube centroids_new(arma::size(centroids));
		arma::vec w(S);
		bool stopit = false;
		int iter = 0;

		while(!stopit && iter < max_iter) {

			// Check interrupt every 100 iter's
			if(iter % 100 == 0) {
				Rcpp::checkUserInterrupt();
			}

			for(int s = 0; s < S; ++s) {
				// Update distance matrix
				for(int k = 0; k < K; ++k) {
					// Integrated weighted distances over n locations
					double dist_sk = 0;
					for(int i = 0; i < n; ++i) {
						dist_sk += std::pow(pdDist_C(M.slice(i * S + s), centroids.slice(k * n + i), metric), (double)2);
					}
					dist_mat(s, k) = dist_weights(s, k) * dist_sk;
				}
				if(m > 1) {
					// Update fuzzy cluster assignments
					if(dist_mat.row(s).min() >= 1E-10) {
						// Update for non-zero distances
						for(int k = 0; k < K; ++k) {
							clust_mat(s, k) = std::pow(dist_mat(s, k), -1 / (m - 1)) /
									arma::accu(arma::pow(dist_mat.row(s), -1 / (m - 1)));
						}
					}
					else {
						// Update for (approximately) zero distances
						clust_mat.row(s).zeros();
						for(int k = 0; k < K; ++k) {
							if(dist_mat(s, k) < 1E-10) {
								clust_mat(s, k) = 1 / arma::accu(dist_mat.row(s) < 1E-10);
							}
						}
					}
				}
			}
			if(m == 1) {
				// Update hard cluster assignments
				arma::uvec min_dist = arma::index_min(dist_mat, 1);
				clust_mat.zeros();
				for(int s = 0; s < S; ++s){
					clust_mat(s, min_dist(s)) = 1;
				}
			}
			// Compute centroids based on new cluster weights
			for(int k = 0; k < K; ++k) {
				// Compute centroid weights
				if(m == 1) {
					w = clust_mat.col(k) / arma::accu(clust_mat.col(k));
				}
				else {
					w = arma::pow(clust_mat.col(k), m) / arma::accu(arma::pow(clust_mat.col(k), m));
				}
				for(int i = 0; i < n; ++i) {
					if(metric == "Riemannian") {
						// Riemannian weighted average
						centroids_new.slice(k * n + i) = pdMean_C_approx(M.slices(i * S, (i + 1) * S - 1), w);
					}
					else {
						// Euclidean weighted average
						arma::cx_mat cent_k = arma::zeros<arma::cx_mat>(d, d);
						for(int s = 0; s < S; ++s){
							cent_k += w(s) * M.slice(i * S + s);
						}
						centroids_new.slice(k * n + i) = cent_k;
					}
				}
			}

			// Revalidate stopping criterion
			double dist_old_new = 0;
			for(arma::uword i = 0; i < centroids.n_slices; ++i) {
				dist_old_new += std::pow(pdDist_C(centroids.slice(i), centroids_new.slice(i), metric), (double)2);
			}
			stopit = (dist_old_new < eps);

			// Update for next iteration
			centroids = centroids_new;
			if(iter >= max_iter) {
				Rcpp::warning("Reached maximum number of iterations in c-medoids algorithm.");
			}
			++iter;
		}

		// Return final cluster assignments and distance matrix
		arma::cube res(S, K, 2);
		res.slice(0) = clust_mat;
		res.slice(1) = dist_mat;
		return res;

		// catch errors
	} catch(std::exception &ex) {
		forward_exception_to_r(ex);
	} catch(...) {
		Rcpp::exception("c++ exception (unknown reason)");
	}
	return arma::zeros<arma::cube>(1,1,1);
}


