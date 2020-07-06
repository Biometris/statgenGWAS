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
  // get maximum number of processes.
  nThr = omp_get_num_procs();
  // restrict to max procs - 1.
  nThr = std::max(nThr - 1, 1);
  // Restrict to thread limit from OMP_THREAD_LIMIT. 
  // Set to 2 by CRAN. INT_MAX if unset.
  nThr = std::min(nThr, omp_get_thread_limit());
  // Restrict to max number of threads from OMP_NUM_THREADS. 
  // Initialized when OpenMP is started.
  nThr = std::min(nThr, omp_get_max_threads());
  if (nCores.isNotNull()) {
    nThr = std::min(IntegerVector(nCores)[0], nThr);
  } 
#endif  
  return nThr;
}
