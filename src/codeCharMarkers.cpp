#include <algorithm>
#include <Rcpp.h>
#include "getThr.h"

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

std::string::size_type countSub(std::string fullString,
                                std::string subString)
{
  int n = 0;
  size_t pos = fullString.find(subString, 0);
  while (pos != std::string::npos)
  {
    n ++;
    pos = fullString.find(subString, pos + 1);
  }
  return n;
}

// [[Rcpp::export]]
NumericMatrix codeCharMarkers(CharacterMatrix& markers,
                              std::vector< std::string > refAlls,
                              bool& minor,
                              double& maxAll,
                              Rcpp::Nullable<Rcpp::IntegerVector> nCores = R_NilValue) 
{
  int nMrk = markers.ncol();
  int nGeno = markers.nrow();
  NumericMatrix out = NumericMatrix( nGeno, nMrk );
  colnames(out) = colnames(markers);
  rownames(out) = rownames(markers);
  std::string refAll;
#ifdef _OPENMP
  int nThr = getThr(nCores);
#pragma omp parallel for private(refAll) num_threads(nThr)
#endif
  for (int i = 0; i < nMrk; i++) 
  {
    refAll = refAlls[i];
    for (int j = 0; j < nGeno; j++) 
    {
      std::string s = Rcpp::as<std::string>(markers(j, i));
      out(j, i) = countSub(s, refAll);
    }
  }
  if (minor) 
  {
    for (int i = 0; i < nMrk; i++) 
    {
      NumericMatrix::Column outCol = out(_, i);
      if (mean(na_omit(outCol)) > maxAll / 2)
      {
        outCol = maxAll - outCol;
      }
    }
  }
  return out;
}


