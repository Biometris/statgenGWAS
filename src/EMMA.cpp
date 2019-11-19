#include <RcppArmadillo.h>

#define _USE_MATH_DEFINES // for C++
#include <cmath>

// Correctly setup the build environment
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' EMMA helper functions
//'
//' Helper functions for computing REML estimates of genetic and residual
//' variance components using the EMMA algorithm.
//'
//' @inheritParams EMMA
//' @param x a q x n covariate matrix, q being the number of covariates and n
//' being the number of genotypes. q has to be at least one (typically an
//' intercept).
//'
//' @noRd
//' @keywords internal
// [[Rcpp::export]]
void emmaEigenR(const arma::mat k,
                const arma::mat x,
                arma::vec &eigVals,
                arma::mat &eigVecs) {
  int n = x.n_rows;
  int q = x.n_cols;
  // Compute n-q non-zero eigenvalues of SHS as defined in eqn. 5 of Kang.
  arma::mat s = arma::eye(n, n);
  if (q == 1) {
    s.for_each( [ n ](mat::elem_type& val) { val -= 1 / double (n); } );
  } else {
    s -= x * solve(x.t() * x, x.t());
  }
  arma::eig_sym(eigVals, eigVecs, s * (k + eye(n, n)) * s);
  eigVals = reverse(eigVals.tail(n - q)) - 1;
  eigVecs = fliplr( eigVecs.tail_cols(n - q) );
}

// [[Rcpp::export]]
arma::vec emmaREMLLL(double logDelta,
                     arma::vec lambda,
                     arma::vec etas1,
                     double n,
                     double t,
                     arma::vec etas2) {
  // Compute the REML LL as in eqn. 7 of Kang.
  double nq = etas1.n_elem + n - t;
  double delta = exp(logDelta);
  lambda += delta;
  arma::vec ll = 0.5 * (nq * (log(nq / (2 * M_PI)) - 1 - log(sum(square(etas1) /
    lambda) + etas2 / delta)) - sum(log(lambda)) + (t - n) * logDelta);
  return ll;
  //return arma::conv_to< std::vector<double>>::from(ll);
}

// [[Rcpp::export]]
double goldenSectionSearch(double upperBound,
                           double center,
                           double lowerBound,
                           double absolutePrecision,
                           arma::vec lambda,
                           arma::vec etas1,
                           double n,
                           double t,
                           arma::vec etas2) {
  double resphi = (3 - sqrt(5)) / 2;
  // a and b are the current bounds; the minimum is between them.
  // c is the center pointer pushed slightly left towards a
  if (abs(upperBound - lowerBound) < absolutePrecision) {
    return (upperBound + lowerBound) / 2;
  }
  //  Create a new possible center, in the area between c and b, pushed against c
  double centerNew = center + resphi * (upperBound - center);
  if (arma::as_scalar(emmaREMLLL(centerNew, lambda, etas1, n, t, etas2)) <
    arma::as_scalar(emmaREMLLL(center, lambda, etas1, n, t, etas2))) {
    return goldenSectionSearch(center, centerNew, upperBound,
                               absolutePrecision, lambda, etas1, n, t, etas2);
  } else {
    return goldenSectionSearch(centerNew, center, lowerBound,
                               absolutePrecision, lambda, etas1, n, t, etas2);
  }
}

// [[Rcpp::export]]
List emmaCPP(arma::vec y,
             arma::mat k,
             arma::mat x,
             int nGrids = 100,
             double uLim = 10,
             double lLim = -10,
             double eps = 1e-3) {
  // Using notation of Kang et al.
  int n = y.n_elem;
  int t = k.n_rows;
  int q = x.n_cols;
  // Define intervals used for computing local optimums.
  double m = nGrids + 1;
  arma::vec logDelta = arma::linspace<vec>(lLim, uLim, nGrids + 1);
  arma::vec delta = exp(logDelta);
  // Compute n-q non-zero eigenvalues and corresponding eigenvectors.
  arma::vec eigVals(t);
  arma::mat eigVecs(t, t);
  emmaEigenR(k, x, eigVals, eigVecs);
  // Define eta as in eqn. 7 of Kang.
  arma::vec etas1 = arma::vectorise(eigVecs.t() * y);
  // Define etas2 to use same optimisation as for Z not NULL.
  arma::vec etas2 = zeros(1);
  // Define lambda as in eqn. 7 of Kang for usage in vectorised form.
  arma::mat lambdas = zeros(n - q , m);
  lambdas.each_col() += eigVals;
  lambdas.each_row() += delta.t();
  arma::mat s1 = 1 / square(lambdas);
  s1.each_col() %= square(etas1);
  arma::mat s2 = 1 / lambdas;
  s2.each_col() %= square(etas1);
  // Compute derivative of LL as in eqn. 9 of Kang for all grid endpoints.
  arma::vec dLL = vectorise(0.5 * delta.t() %
    ((n - q) * sum(s1) / sum(s2) - sum(1 / lambdas)));
  // Find optimum of LL
  arma::vec optLogDelta(m);
  arma::vec optLL(m);
  int rel = 0;
  // Check first item in dLL. If < eps include LL value as possible optimum.
  if (dLL(0) < eps) {
    optLogDelta(rel) = lLim;
    optLL(rel) = arma::as_scalar(emmaREMLLL(lLim, eigVals, etas1, 0, 0, etas2));
    rel ++;
  }
  // Check last item in dLL. If > - eps include LL value as possible optimum.
  if (dLL(m - 1) > -eps) {
    optLogDelta(rel) = uLim;
    optLL(rel) = arma::as_scalar(emmaREMLLL(uLim, eigVals, etas1, 0, 0, etas2));
    rel ++;
  }
  // If derivative changes sign on an interval, compute local optimum for LL
  // on that interval and add it to possible optimums.
  for (unsigned int i = 0; i < m - 2; i++) {
    if ((dLL(i) > 0 && dLL(i + 1) < 0) || dLL(i) * dLL(i + 1) < pow(eps, 2)) {
      optLogDelta(rel) = goldenSectionSearch(logDelta(i + 1),
                  (-1 + sqrt(5)) / 2 * (logDelta(i) + logDelta(i + 1)),
                  logDelta(i), eps, eigVals, etas1, 0, 0, etas2);
      optLL(rel) = arma::as_scalar(emmaREMLLL(optLogDelta(rel), eigVals, etas1,
                                   0, 0, etas2));
      rel ++;
    }
  }
  optLL = optLL.head(rel);
  optLogDelta = optLogDelta.head(rel);
  // Compute absolute LL maximum from possible optimums.
  double maxDelta = exp(optLogDelta( optLL.index_max() ));
  // Compute variance components.
  double maxVg = sum(square(etas1) / (eigVals + maxDelta)) / (n - q);
  double maxVe = maxVg * maxDelta;
  arma::mat vcovMatrix = maxVg * k + maxVe * eye(t, t);
  return List::create(_["maxVg"] = maxVg,
                      _["maxVe"] = maxVe,
                      _["vcovMatrix"] = vcovMatrix);
}
