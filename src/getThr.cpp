#include <Rcpp.h>
using namespace Rcpp;

// Add a flag to enable OpenMP at compile time
// [[Rcpp::plugins(openmp)]]

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
int getThr(Rcpp::Nullable<Rcpp::IntegerVector> nCores = R_NilValue) {
  int nThr = 1;
#ifdef _OPENMP
  int maxThr = omp_get_max_threads();
  if (nCores.isNotNull()) {
    nThr = std::min(IntegerVector(nCores)[0], maxThr);
  } else {
    nThr = maxThr - 1;
  }
#endif
  return nThr;
}
