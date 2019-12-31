#include <RcppArmadillo.h>

// Correctly setup the build environment
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Reduce the kinship matrix
//'
//' The kinship matrix is reduced using nPca eigenvectors of K.
//'
//' @inheritParams runSingleTraitGwas
//'
//' @param nPca An integer, the number of eigenvectors used for reducing the
//' kinship matrix.
//'
//' @return The reduced kinship matrix
//'
//' @noRd
//' @keywords internal
// [[Rcpp::export]]
arma::mat reduceKinship(const arma::mat K,
                        const int nPca) {
  arma::vec eigVals( K.n_cols );
  arma::mat eigVecs( size(K) );
  arma::eig_sym(eigVals, eigVecs, K);
  // fliplr and reverse are needed because eigenVals and eigenVecs are returned
  // sorted in ascending order.
  eigVecs = fliplr( eigVecs.tail_cols(nPca) );
  arma::mat S = arma::diagmat( reverse(eigVals.tail(nPca)) );
  return eigVecs * S * eigVecs.t();
}

// [[Rcpp::export]]
arma::mat nearestPD(arma::mat x,
                    const bool corr = false,
                    const bool keepDiag = false,
                    const bool do2eigen = true,
                    const bool doSym = false,
                    const bool doDykstra = true,
                    const double eigTol = 1e-6,
                    const double convTol = 1e-7,
                    const double posdTol = 1e-8,
                    const int maxIter = 100) {
  if (!x.is_symmetric()) {
    x = (x + x.t()) / 2;
  }
  int n = x.n_cols;
  arma::vec diagX0 = zeros(n);
  arma::mat DykS = zeros( size(x) );
  if (keepDiag) {
    diagX0 = x.diag();
  }
  arma::mat xUpdate = x;
  int iter = 0;
  bool converged = false;
  double conv = std::numeric_limits<double>::max();
  arma::mat r = zeros( size(x) );
  while (iter < maxIter && !converged) {
    arma::mat y = xUpdate;
    arma::vec eigVals(n);
    arma::mat eigVecs( size(x) );
    if (doDykstra) {
      r = y - DykS;
      arma::eig_sym(eigVals, eigVecs, r);
    } else {
      arma::eig_sym(eigVals, eigVecs, y);
    }
    uvec p = find (eigVals > eigTol * as_scalar( eigVals.tail(1)) );
    // if (!any(p))
    //   stop("Matrix seems negative semi-definite")
    eigVecs = eigVecs.cols( p );
    arma::mat eigVecsT = eigVecs.t();
    eigVecs.each_row() %= eigVals( p ).t();
    xUpdate = eigVecs * eigVecsT;
    if (doDykstra) {
      DykS = xUpdate - r;
    }
    if (doSym) {
      xUpdate = (xUpdate + xUpdate.t()) / 2;
    }
    if (corr) {
      xUpdate.diag().ones();
    } else if (keepDiag) {
      xUpdate.diag() = diagX0;
    }
    conv = norm(y - xUpdate, "inf") / norm(y, "inf");
    iter ++;
    converged = (conv <= convTol);
  } // end while
  if (!converged) {
     warning("'nearestPD()' did not converge in " + std::to_string(iter) +
       " iterations.\n");
  }
  if (do2eigen) {
    arma::vec eigVals(n);
    arma::mat eigVecs( size(x) );
    arma::eig_sym(eigVals, eigVecs, xUpdate);
    double eps = posdTol * fabs(as_scalar(eigVals.tail(1)));
    if (eigVals(0) < eps) {
      eigVals = clamp(eigVals, eps, eigVals.max());
      arma::vec oDiag = xUpdate.diag();
      arma::mat eigVecsT = eigVecs.t();
      eigVecsT.each_col() %= eigVals;
      xUpdate = eigVecs * eigVecsT;
      arma::vec dVec = sqrt(clamp(oDiag, eps, oDiag.max()) / xUpdate.diag());
      xUpdate.each_row() %= dVec.t();
      xUpdate.each_col() %= dVec;
    }
    if (corr) {
      xUpdate.diag().ones();
    } else if (keepDiag) {
      xUpdate.diag() = diagX0;
    }
  }
  return xUpdate;
}



