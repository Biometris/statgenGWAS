#ifndef GETTHR_H_
#define GETTHR_H_

#include "Rcpp.h"

int getThr(Rcpp::Nullable<Rcpp::IntegerVector> nCores = R_NilValue);

#endif
