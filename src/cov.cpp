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

//printf("n=%d, Ntau=%d, Nsigma=%d, Nphi=%d\n", n, Ntau, Nsigma, Nphi);

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
						(phi[R[i]]*phi[R[j]]/rt_dets(R[i], R[j])) * 
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

// [[Rcpp::export]]
NumericMatrix calc_ns_cov_2phi_Rcpp(
	NumericVector tau, NumericVector sigma, NumericVector phi1, NumericVector phi2,
	int Nr, NumericVector R, NumericMatrix S
) {

	int n      = S.nrow();      // number of locations
	int Ntau   = tau.size();    // number of nugget params
	int Nsigma = sigma.size();  // number of partial sill params
	int Nphi1  = phi1.size();   // number of range params
	int Nphi2  = phi2.size();

printf("n=%d (%d), Ntau=%d, Nsigma=%d, Nphi1=%d, Nphi2=%d\n", n, S.ncol(), Ntau, Nsigma, Nphi1, Nphi2);

	// fill in covariance matrix
	NumericMatrix Sigma(n, n);

	int i,j;
	double d1,d2;

	// start with range
	if (Nphi1 == 1 && Nphi2 == 1) {
		// range matrix doesn't vary
		double sqr_phi1 = pow(phi1[0],2.0);
		double sqr_phi2 = pow(phi2[0],2.0);

		for (i = 0; i < n; i++) for (j = i; j < n; j++) {
			d1 = pow(S(i,0)-S(j,0), 2.0);
			d2 = pow(S(i,1)-S(j,1), 2.0);

			Sigma(i,j) = exp( -sqrt( d1/sqr_phi1 + d2/sqr_phi2 ) );
		}
	} else {
		// varying range

/*
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
*/

		NumericVector phi_prods(Nr);
		NumericVector phi1_sqrs(Nr);
		NumericVector phi2_sqrs(Nr);
		for (i = 0; i < Nr; i++) {
			phi_prods[i] = sqrt(phi1[i]*phi2[i]);
			phi1_sqrs[i] = pow(phi1[i],2);
			phi2_sqrs[i] = pow(phi2[i],2);
		}

		for (i = 0; i < n; i++) {
			for (j = i; j < n; j++) {
				d1 = pow(S(i,0)-S(j,0), 2.0);
				d2 = pow(S(i,1)-S(j,1), 2.0);
//if (i<10&j<10) printf("[%d,%d] %.2f, %.2f (%.2f, %.2f), (%.2f, %.2f)\n", i+1,j+1, d1, d2, S(i,0), S(i,1), S(j,0), S(j,1));

				if (R[i] == R[j]) {
					// same region
					Sigma(i,j) = exp( -sqrt( d1/pow(phi1[R[i]],2.0) + d2/pow(phi2[R[i]],2.0) ) );
				} else {
					// different regions
					Sigma(i,j) =
						(2.0*phi_prods[R[i]] * phi_prods[R[j]] / sqrt( (phi1_sqrs[R[i]]+phi1_sqrs[R[j]])*(phi2_sqrs[R[i]]+phi2_sqrs[R[j]]) )) *
						exp( -sqrt( 2.0*(d1/(phi1_sqrs[R[i]]+phi1_sqrs[R[j]]) + d2/(phi2_sqrs[R[i]]+phi2_sqrs[R[j]])) ) );

/*
						2*sqrt(phi1[R[i]]*phi2[R[i]] * phi1[R[j]]*phi2[R[j]] /
						( (pow(phi1[R[i]],2)+pow(phi1[R[j]],2)) * (pow(phi2[R[i]],2)+pow(phi2[R[j]],2)) ) *
						exp( -sqrt( d1/(pow(phi1[R[i]],2)+pow(phi1[R[j]],2)) + d2/(pow(phi2[R[i]],2)+pow(phi2[R[j]],2)) ) );
*/
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

// [[Rcpp::export]]
NumericMatrix calc_ns_cov_angle_Rcpp(
	NumericVector tau, NumericVector sigma, NumericVector phi1, NumericVector phi2, NumericVector rho,
	int Nr, NumericVector R, NumericMatrix S
) {

	int n      = S.nrow();      // number of locations
	int Ntau   = tau.size();    // number of nugget params
	int Nsigma = sigma.size();  // number of partial sill params
	int Nphi1  = phi1.size();   // number of range params
	int Nphi2  = phi2.size();
	int Nrho   = rho.size();

printf("n=%d (%d), Ntau=%d, Nsigma=%d, Nphi1=%d, Nphi2=%d, Nrho=%d\n", n, S.ncol(), Ntau, Nsigma, Nphi1, Nphi2, Nrho);

	// fill in covariance matrix
	NumericMatrix Sigma(n, n);

	int i,j;
	double d1,d2,d1_sqr,d2_sqr;
	double detOmega;
	NumericMatrix invOmega(2, 2);

	// start with range
	if (Nphi1 == 1 && Nphi2 == 1 && Nrho == 1) {
		// range matrix doesn't vary
		double sqr_phi1 = pow(phi1[0],2.0);
		double sqr_phi2 = pow(phi2[0],2.0);

		// invert Omega
		detOmega = sqr_phi1*sqr_phi2 - pow(rho[0]*phi1[0]*phi2[0], 2.0);
		invOmega(0,0) = sqr_phi1/detOmega;
		invOmega(1,1) = sqr_phi2/detOmega;
		invOmega(0,1) = -rho[0]*phi1[0]*phi2[0]/detOmega;

		for (i = 0; i < n; i++) for (j = i; j < n; j++) {
			d1 = S(i,0)-S(j,0);
			d2 = S(i,1)-S(j,1);
			d1_sqr = pow(d1, 2.0);
			d2_sqr = pow(d2, 2.0);

			Sigma(i,j) = exp( -sqrt( d1_sqr*invOmega(0,0) + d2_sqr*invOmega(1,1) + 2*d1*d2*invOmega(0,1) ) );
		}
	} else {
		// varying range

/*
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
*/

		NumericVector phi_prods(Nr);
		NumericVector rt_phi_prods(Nr);
		NumericVector phi1_sqrs(Nr);
		NumericVector phi2_sqrs(Nr);
		NumericVector all_prods1(Nr);
		NumericVector all_prods2(Nr);
		for (i = 0; i < Nr; i++) {
			phi_prods[i]    = phi1[i]*phi2[i];
			rt_phi_prods[i] = sqrt(phi_prods[i]);
			phi1_sqrs[i]    = pow(phi1[i],2.0);
			phi2_sqrs[i]    = pow(phi2[i],2.0);
			all_prods1[i]   = phi1[i]*phi2[i]*sqrt(1.0-pow(rho[i],2.0));
			all_prods2[i]   = rho[i]*phi1[i]*phi2[i];
		}

		double rt2 = sqrt(2.0);

		for (i = 0; i < n; i++) {
			for (j = i; j < n; j++) {
				d1 = S(i,0)-S(j,0);
				d2 = S(i,1)-S(j,1);
				d1_sqr = pow(d1, 2.0);
				d2_sqr = pow(d2, 2.0);
//if (i<10&j<10) printf("[%d,%d] %.2f, %.2f (%.2f, %.2f), (%.2f, %.2f)\n", i+1,j+1, d1, d2, S(i,0), S(i,1), S(j,0), S(j,1));

				// invert Omega
				detOmega = (phi1_sqrs[R[i]]+phi1_sqrs[R[j]])*(phi2_sqrs[R[i]]+phi2_sqrs[R[j]]) - pow(all_prods2[R[i]]+all_prods2[R[j]], 2.0);
				invOmega(0,0) = (phi1_sqrs[R[i]]+phi1_sqrs[R[j]])/detOmega;
				invOmega(1,1) = (phi2_sqrs[R[i]]+phi2_sqrs[R[j]])/detOmega;
				invOmega(0,1) = -(all_prods2[R[i]]+all_prods2[R[j]])/detOmega;

				if (R[i] == R[j]) {
					// same region
					Sigma(i,j) = exp( -rt2*sqrt( d1_sqr*invOmega(0,0) + d2_sqr*invOmega(1,1) + 2.0*d1*d2*invOmega(0,1) ) );
				} else {
					// different regions
					Sigma(i,j) =
						(2.0 * sqrt( all_prods1[R[i]]*all_prods1[R[j]]/( (phi1_sqrs[R[i]]+phi1_sqrs[R[j]])*(phi2_sqrs[R[i]]+phi2_sqrs[R[j]]) - pow(all_prods2[R[i]]+all_prods2[R[j]], 2.0) ) ) ) *
						exp(-rt2*sqrt( d1_sqr*invOmega(0,0) + d2_sqr*invOmega(1,1) + 2.0*d1*d2*invOmega(0,1) ));

/*
double mid = d1_sqr*invOmega(0,0) + d2_sqr*invOmega(1,1) + 2.0*d1*d2*invOmega(0,1);
if (mid < 0) { printf("Doh: (%d,%d)=%.3f: d1^2=%.3f, d2^2=%.3f, d1=%.3f, d2=%.3f, d1*d2=%.3f\n", mid, i, j, d1_sqr, d2_sqr, d1, d2, d1*d2); break; }
*/
/*
						2*sqrt(phi1[R[i]]*phi2[R[i]] * phi1[R[j]]*phi2[R[j]] /
						( (pow(phi1[R[i]],2)+pow(phi1[R[j]],2)) * (pow(phi2[R[i]],2)+pow(phi2[R[j]],2)) ) *
						exp( -sqrt( d1/(pow(phi1[R[i]],2)+pow(phi1[R[j]],2)) + d2/(pow(phi2[R[i]],2)+pow(phi2[R[j]],2)) ) );
*/
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
