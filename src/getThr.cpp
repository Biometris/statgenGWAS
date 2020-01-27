#include <Rcpp.h>
using namespace Rcpp;

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::export]]
int getThr(Rcpp::Nullable<Rcpp::IntegerVector> nCores = R_NilValue) {
  int nThr = 1;
#ifdef _OPENMP
  // get thread limit from OMP_THREAD_LIMIT. Set to 2 by CRAN. INT_MAX if unset.
  int thrLim = omp_get_thread_limit();  
  // get max number of threads from OMP_NUM_THREADS. 
  // Initialized when OpenMP is started.
  int maxThr = omp_get_max_threads();
  int maxThrUsed = std::min(thrLim, maxThr);
  if (nCores.isNotNull()) {
    nThr = std::min(IntegerVector(nCores)[0], maxThrUsed);
  } else {
    nThr = maxThrUsed - 1;
  }
#endif
  return nThr;
}
