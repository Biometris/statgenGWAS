#include <RcppArmadillo.h>
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
arma::mat astleCPP(arma::mat x,
                   Rcpp::Nullable<Rcpp::NumericVector> denom = R_NilValue) {
  // Remove markers with variance 0.
  x = x.cols( find (var(x) > 0) );
  // Scale X.
  arma::rowvec p = sum(x) / (2 * x.n_rows);
  x.each_row() -= 2 * p;
  x.each_row() /= sqrt(2 * p % (1 - p));
  // Compute denominator.
  double denominator;
  if (denom.isNull()) {
    denominator = x.n_cols;
  } else {
    denominator = Rcpp::as<double>(denom);
  }
  return x * x.t() / denominator;
}

// [[Rcpp::export]]
arma::mat IBSCPP(arma::mat x,
                 Rcpp::Nullable<Rcpp::NumericVector> denom = R_NilValue) {
  // Remove markers with variance 0.
  x = x.cols( find (var(x) > 0) );
  // Compute denominator.
  double denominator;
  if (denom.isNull()) {
    denominator = x.n_cols;
  } else {
    denominator = Rcpp::as<double>(denom);
  }
  return (x * x.t() + (1 - x) * (1 - x).t()) / denominator;
}

// [[Rcpp::export]]
arma::mat vanRadenCPP(arma::mat x,
                      Rcpp::Nullable<Rcpp::NumericVector> denom = R_NilValue) {
  // Remove markers with variance 0.
  x = x.cols( find (var(x) > 0) );
  // Scale X.
  arma::rowvec p = sum(x) / (2 * x.n_rows);
  x.each_row() -= 2 * p;
  // Compute denominator.
  double denominator;
  if (denom.isNull()) {
    denominator = 2 * sum(p % (1 - p));
  } else {
    denominator = Rcpp::as<double>(denom);
  }
  return x * x.t() / denominator;
}
