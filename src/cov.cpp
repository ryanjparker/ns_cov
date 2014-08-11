#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double meanC(NumericVector x) {
  int n = x.size();
  double total = 0;

  for(int i = 0; i < n; ++i) {
    total += x[i];
  }
  return total / n;
}

// [[Rcpp::export]]
double sum_diag_mm_Rcpp(NumericMatrix A, NumericMatrix B) {
	int n = A.nrow();
	int i,j;
	double s = 0;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			s += A(i,j) * B(j,i);
		}
	}

	return(s);
}

// [[Rcpp::export]]
NumericMatrix calc_ns_cov_Rcpp(
	NumericVector tau, NumericVector sigma, NumericVector phi,
	int Nr, NumericVector R, NumericMatrix D2
) {

	int n      = D2.nrow();     // number of locations
	int Ntau   = tau.size();    // number of nugget params
	int Nsigma = sigma.size();  // number of partial sill params
	int Nphi   = phi.size();    // number of range params

	// fill in covariance matrix
	NumericMatrix Sigma(n, n);

	int i,j;

	// start with phi
	if (Nphi == 1) {
		// single range
		for (i = 0; i < n; i++) for (j = i; j < n; j++) Sigma(i,j) = exp(-sqrt(D2(i,j))/phi[0]);
	} else {
		// varying range

		// compute root determinants and inverses between regions
		NumericMatrix rt_dets(Nr, Nr);
		NumericMatrix invs(Nr, Nr);
		double S2;
		for (i = 0; i < Nr; i++) {
			for (j = i; j < Nr; j++) {
				S2 = 0.5 * (pow(phi[i], 2) + pow(phi[j], 2));
				rt_dets(i, j) = rt_dets(j, i) = S2;
				invs(i, j) = invs(j, i) = 1/S2;
			}
		}

		for (i = 0; i < n; i++) {
			for (j = i; j < n; j++) {
				if (R[i] == R[j]) {
					// same region
					Sigma(i,j) = exp(-sqrt(D2(i,j))/phi[R[i]]);
				} else {
					// different regions
					Sigma(i,j) =
						phi[R[i]]*phi[R[j]]/rt_dets(R[i], R[j]) * 
						exp(-sqrt( D2(i,j) * invs(R[i], R[j]) ));
				}
			}
		}
	}

	// scale by sigma
	if (Nsigma == 1) {
		// single partial sill
		for (i = 0; i < n; i++) for (j = i; j < n; j++) Sigma(i,j) *= pow(sigma[0], 2);
	} else {
		// varying partial sill
		for (i = 0; i < n; i++) for (j = i; j < n; j++) Sigma(i,j) *= sigma[R[i]]*sigma[R[j]];
	}

	// add nugget
	if (Ntau == 1) {
		// single nugget
		for (i = 0; i < n; i++) Sigma(i,i) += tau[0];
	} else {
		// varying nugget
		for (i = 0; i < n; i++) Sigma(i,i) += tau[R[i]];
	}

	// fill in lower triangle
	for (i = 0; i < n; i++) for (j = i+1; j < n; j++) Sigma(j,i) = Sigma(i,j);

	return(Sigma);
}

/*

// [[Rcpp::export]]
NumericMatrix ce_cov_Rcpp(NumericVector theta, NumericMatrix X) {
	int n      = X.nrow();
	int Ntheta = theta.size();

	NumericMatrix Sigma(n, n);

	int i,j,k;

	for (i = 0; i < n; i++) {
		for (j = i; j < n; j++) {
			Sigma(i,j) = 0;
			for (k = 0; k < Ntheta; k++) {
				Sigma(i,j) += theta[k]*pow(X(i,k)-X(j,k), 2);
			}
			Sigma(i,j) = exp(-Sigma(i,j));
			Sigma(j,i) = Sigma(i,j);
		}
	}

	return(Sigma);
}

// [[Rcpp::export]]
NumericMatrix ce_partial_Rcpp(int e, NumericVector theta, NumericMatrix Sigma, NumericMatrix X) {
	int n      = X.nrow();
	int Ntheta = theta.size();

	NumericMatrix P(n, n);

	int i,j,k;
	k = e-1;

	double et = exp(theta[k]);

	for (i = 0; i < n; i++) {
		for (j = i; j < n; j++) {
			P(i,j) = -pow(X(i,k)-X(j,k), 2) * et * Sigma(i,j);
			P(j,i) = P(i,j);
		}
	}

	return(P);
}

// [[Rcpp::export]]
NumericMatrix ce_full_pred_X_cov(NumericMatrix X, NumericMatrix Xobs, NumericVector theta) {
	int Npred  = X.nrow();
	int Nobs   = Xobs.nrow();
	int Ntheta = theta.size();

	NumericMatrix Sigma(Npred, Nobs);

	int i,j,k;

	for (i = 0; i < Npred; i++) {
		for (j = 0; j < Nobs; j++) {
			Sigma(i,j) = 0;
			for (k = 0; k < Ntheta; k++) {
				Sigma(i,j) += theta[k]*pow(X(i,k)-Xobs(j,k), 2);
			}
			Sigma(i,j) = exp(-Sigma(i,j));
		}
	}

	return(Sigma);
}

*/
