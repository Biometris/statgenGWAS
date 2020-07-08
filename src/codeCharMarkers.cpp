#include <algorithm>
#include <Rcpp.h>
#include "getThr.h"

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// Helper function for counting the number of times a sub string appears in a 
// full string. 
// Used for counting the number of alleles. 
int countSub(std::string fullString,
             std::string subString) {
  // Initialize n to 0, results in returning 0 if no sub string found.
  int n = 0;
  // Find first position of sub string.
  size_t pos = fullString.find(subString, 0);
  // While not at the end of the string continue looking for the next occurrence
  // of sub string and increase n in the process.
  while (pos != std::string::npos) {
    n ++;
    pos = fullString.find(subString, pos + 1);
  }
  return n;
}

// Function for coding marker matrix with character entries.

// [[Rcpp::export]]
List codeCharMarkers(CharacterMatrix& markers,
                     std::vector< std::string > refAlls,
                     bool& minor,
                     Rcpp::Nullable<Rcpp::IntegerVector> nCores = R_NilValue) {
  // Initialize output matrix;
  int nMrk = markers.ncol();
  int nGeno = markers.nrow();
  NumericMatrix out = NumericMatrix( nGeno, nMrk );
  colnames(out) = colnames(markers);
  rownames(out) = rownames(markers);
  // Determine the number of times punctuation occurs in the first non NA 
  // entry in the first marker column.
  // There always is such an entry since column containing only NA have been
  // removed earlier.
  std::string firstAll0 = Rcpp::as<std::string>(na_omit(markers.column(0))[0]);
  const int nPunct = std::count_if(firstAll0.begin(), firstAll0.end(), ::ispunct);
  // Determine the allele length from the total length and the number of times
  // punctuation is found. 
  // This assumes an input where all alleles have the same length. Other inputs
  // don't make sense, so are just ignored.
  const int allLength = (firstAll0.length() - nPunct) / (nPunct + 1);
  // Compute the maximum number of alleles that might occur.
  double maxAll;
  if (nPunct > 0) {
    // with punctuation, maximum length is one more than number of punctuations.
    maxAll = nPunct + 1;
  } else {
    // no punctuation, maximum length is the full length.
    maxAll = allLength;
  }
#ifdef _OPENMP
  int nThr = getThr(nCores);
#pragma omp parallel for num_threads(nThr)
#endif
  for (int i = 0; i < nMrk; i++) {
    std::string refAll;
    if (minor) {
      // Get the first non missing entry for the current marker.
      // Using na_omit here crashes when using openMP so working around that.
      int j = 0;
      while (markers(j, i) == NA_STRING) {
        j++;
      }  
      std::string firstAll = Rcpp::as<std::string>(markers(j, i));
      // Extract the reference allele.
      if (nPunct > 0) {
        // With punctuation, use computed allele length.
        refAll = firstAll.substr(0, allLength);
      } else {
        // No punctuation, the first character.
        refAll = firstAll.at(0);
      }
    } else {
      // Not using the minor allele as reference allele.
      // Get the reference allele from the input.
      refAll = refAlls[i];
    }
    // Initialize column sum and number for computing average.
    double colSum = 0, colNum = 0;
    for (int j = 0; j < nGeno; j++) {
      if (!Rcpp::CharacterVector::is_na(markers(j, i))) {
        // Get the current marker entry.
        std::string s = Rcpp::as<std::string>(markers(j, i));
        // Count the number of times the reference allele occurs in the entry.
        int nAll = countSub(s, refAll);
        // Add to the output and increase column sum and number.
        out(j, i) = nAll;
        colSum = colSum + nAll;
        colNum++;
      } else {
        // Missing entry. Set output to missing as well.
        out(j, i) = NA_REAL;
      }
    }
    // When using minor allele as reference allele, column mean can never be 
    // higher than half of maximum number of alleles.
    // If this is the case recompute values by subtracting them from maximum
    // number of alleles.
    if (minor && (colSum / colNum) > (maxAll / 2)) {
      for (int j = 0; j < nGeno; j++) {
        if (!Rcpp::NumericVector::is_na(out(j, i))) {
          out(j, i) = maxAll - out(j, i);
        }
      }
    }
  }
  // Create output.
  // maxAll is needed for other steps in the process.
  List res = List::create(Named("markersRecoded") = out, 
                          Named("maxAll") = maxAll);
  return res;
}
