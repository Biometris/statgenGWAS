#include <RcppArmadillo.h>
#include "getThr.h"

// Correctly setup the build environment
// [[Rcpp::depends(RcppArmadillo)]]

// Add a flag to enable OpenMP at compile time
// [[Rcpp::plugins(openmp)]]

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List fastGLSCPP(const arma::mat &X,
                const arma::vec &y,
                const arma::mat &sigma,
                Rcpp::Nullable<Rcpp::NumericVector> size_param = R_NilValue,
                Rcpp::Nullable<Rcpp::IntegerVector> nCores = R_NilValue) {
  // Get number of genotypes and markers.
  unsigned int p = X.n_cols;
  unsigned int n = X.n_rows;
  // Define covs as intercept.
  arma::mat covs = ones<mat>(n, 1);
  if (size_param.isNotNull()) {
    // Add other covariates.
    covs.insert_cols(1, Rcpp::as<arma::mat>(size_param));
  }
  // Get number of covariates (including intercept).
  unsigned int nCov = covs.n_cols;
  // Compute inverse of t(M).
  arma::mat mt = chol(sigma, "lower").i();
  // Pre-multiply the snp-matrix with t(M).
  arma::mat tMX = mt * X;
  // Pre-multiply the phenotype (y) with t(M).
  arma::mat tMy = mt * y;
  // pre-multiply the intercept and covariates with t(M).
  arma::mat tMfixCovs = mt * covs;
  // Compute residuals and RSS over all markers.
  arma::mat Q, R;
  arma::qr_econ(Q, R, tMfixCovs);
  arma::vec ResEnv = tMy - Q * Q.t() * tMy;
  double RSSEnv = accu(square(ResEnv));
  // Matrix cookbook, 3.2.6 Rank-1 update of inverse of inner product.
  arma::mat A = inv_sympd(tMfixCovs.t() * tMfixCovs);
  arma::rowvec vv = sum(square(tMX));
  arma::mat vX = tMfixCovs.t() * tMX;
  arma::colvec nn = 1 / (vv - sum(vX % (A * vX))).t();
  arma::mat tvXA = vX.t() * A;
  tvXA.each_row() %= (tMfixCovs.t() * tMy).t();
  arma::colvec betaVec = nn % (tMX.t() * tMy - sum(tvXA, 1));
  // QR decomposition of covariates.
  arma::mat tMQtQ = (mt.t() * (eye<mat>(n, n) - Q * Q.t())).t();
  NumericVector pVal = NumericVector(p);
  NumericVector RLR2 = NumericVector(p);
  // Compute RSS per marker
  arma::mat tX = tMQtQ * X;
  int nThr = getThr(nCores);
#pragma omp parallel for num_threads(nThr)
  for (unsigned int i = 0; i < tX.n_cols; i++) {
    arma::mat Qx, Rx;
    arma::qr_econ(Qx, Rx, tX.col(i));
    double RSSx = accu(square(ResEnv - Qx * Qx.t() * ResEnv));
    double fValx = (RSSEnv - RSSx) / RSSx * (n - 1 - nCov);
    pVal(i) = R::pf(fValx, 1.0, n - 1 - nCov, false, false);
    RLR2(i) = 1 - exp((RSSx - RSSEnv) / n);
  }
  return List::create(_["beta"] = betaVec,
                      _["betaSe"] = sqrt(nn),
                      _["pVal"] = pVal,
                      _["RLR2"] = RLR2);
}


