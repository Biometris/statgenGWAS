// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// emmaEigenR
void emmaEigenR(const arma::mat k, const arma::mat x, arma::vec& eigVals, arma::mat& eigVecs);
RcppExport SEXP _statgenGWAS_emmaEigenR(SEXP kSEXP, SEXP xSEXP, SEXP eigValsSEXP, SEXP eigVecsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type k(kSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type eigVals(eigValsSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type eigVecs(eigVecsSEXP);
    emmaEigenR(k, x, eigVals, eigVecs);
    return R_NilValue;
END_RCPP
}
// emmaREMLLL
arma::vec emmaREMLLL(double logDelta, arma::vec lambda, arma::vec etas1, double n, double t, arma::vec etas2);
RcppExport SEXP _statgenGWAS_emmaREMLLL(SEXP logDeltaSEXP, SEXP lambdaSEXP, SEXP etas1SEXP, SEXP nSEXP, SEXP tSEXP, SEXP etas2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type logDelta(logDeltaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type etas1(etas1SEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type etas2(etas2SEXP);
    rcpp_result_gen = Rcpp::wrap(emmaREMLLL(logDelta, lambda, etas1, n, t, etas2));
    return rcpp_result_gen;
END_RCPP
}
// goldenSectionSearch
double goldenSectionSearch(double upperBound, double center, double lowerBound, double absolutePrecision, arma::vec lambda, arma::vec etas1, double n, double t, arma::vec etas2);
RcppExport SEXP _statgenGWAS_goldenSectionSearch(SEXP upperBoundSEXP, SEXP centerSEXP, SEXP lowerBoundSEXP, SEXP absolutePrecisionSEXP, SEXP lambdaSEXP, SEXP etas1SEXP, SEXP nSEXP, SEXP tSEXP, SEXP etas2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type upperBound(upperBoundSEXP);
    Rcpp::traits::input_parameter< double >::type center(centerSEXP);
    Rcpp::traits::input_parameter< double >::type lowerBound(lowerBoundSEXP);
    Rcpp::traits::input_parameter< double >::type absolutePrecision(absolutePrecisionSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type etas1(etas1SEXP);
    Rcpp::traits::input_parameter< double >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type etas2(etas2SEXP);
    rcpp_result_gen = Rcpp::wrap(goldenSectionSearch(upperBound, center, lowerBound, absolutePrecision, lambda, etas1, n, t, etas2));
    return rcpp_result_gen;
END_RCPP
}
// emmaCPP
List emmaCPP(arma::vec y, arma::mat k, arma::mat x, int nGrids, double uLim, double lLim, double eps);
RcppExport SEXP _statgenGWAS_emmaCPP(SEXP ySEXP, SEXP kSEXP, SEXP xSEXP, SEXP nGridsSEXP, SEXP uLimSEXP, SEXP lLimSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type nGrids(nGridsSEXP);
    Rcpp::traits::input_parameter< double >::type uLim(uLimSEXP);
    Rcpp::traits::input_parameter< double >::type lLim(lLimSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(emmaCPP(y, k, x, nGrids, uLim, lLim, eps));
    return rcpp_result_gen;
END_RCPP
}
// fastGLSCPP
arma::mat fastGLSCPP(const arma::mat& X, const arma::vec& y, const arma::mat& sigma, Rcpp::Nullable<Rcpp::NumericVector> size_param, Rcpp::Nullable<Rcpp::IntegerVector> nCores);
RcppExport SEXP _statgenGWAS_fastGLSCPP(SEXP XSEXP, SEXP ySEXP, SEXP sigmaSEXP, SEXP size_paramSEXP, SEXP nCoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type size_param(size_paramSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::IntegerVector> >::type nCores(nCoresSEXP);
    rcpp_result_gen = Rcpp::wrap(fastGLSCPP(X, y, sigma, size_param, nCores));
    return rcpp_result_gen;
END_RCPP
}
// getThr
int getThr(Rcpp::Nullable<Rcpp::IntegerVector> nCores);
RcppExport SEXP _statgenGWAS_getThr(SEXP nCoresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::IntegerVector> >::type nCores(nCoresSEXP);
    rcpp_result_gen = Rcpp::wrap(getThr(nCores));
    return rcpp_result_gen;
END_RCPP
}
// astleCPP
arma::mat astleCPP(arma::mat x, Rcpp::Nullable<Rcpp::NumericVector> denom);
RcppExport SEXP _statgenGWAS_astleCPP(SEXP xSEXP, SEXP denomSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type denom(denomSEXP);
    rcpp_result_gen = Rcpp::wrap(astleCPP(x, denom));
    return rcpp_result_gen;
END_RCPP
}
// IBSCPP
arma::mat IBSCPP(arma::mat x, Rcpp::Nullable<Rcpp::NumericVector> denom);
RcppExport SEXP _statgenGWAS_IBSCPP(SEXP xSEXP, SEXP denomSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type denom(denomSEXP);
    rcpp_result_gen = Rcpp::wrap(IBSCPP(x, denom));
    return rcpp_result_gen;
END_RCPP
}
// vanRadenCPP
arma::mat vanRadenCPP(arma::mat x, Rcpp::Nullable<Rcpp::NumericVector> denom);
RcppExport SEXP _statgenGWAS_vanRadenCPP(SEXP xSEXP, SEXP denomSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type denom(denomSEXP);
    rcpp_result_gen = Rcpp::wrap(vanRadenCPP(x, denom));
    return rcpp_result_gen;
END_RCPP
}
// matrixRoot
arma::mat matrixRoot(const arma::mat x);
RcppExport SEXP _statgenGWAS_matrixRoot(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(matrixRoot(x));
    return rcpp_result_gen;
END_RCPP
}
// reduceKinship
arma::mat reduceKinship(const arma::mat K, const int nPca);
RcppExport SEXP _statgenGWAS_reduceKinship(SEXP KSEXP, SEXP nPcaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type K(KSEXP);
    Rcpp::traits::input_parameter< const int >::type nPca(nPcaSEXP);
    rcpp_result_gen = Rcpp::wrap(reduceKinship(K, nPca));
    return rcpp_result_gen;
END_RCPP
}
// nearestPD
arma::mat nearestPD(arma::mat x, const bool corr, const bool keepDiag, const bool do2eigen, const bool doSym, const bool doDykstra, const double eigTol, const double convTol, const double posdTol, const int maxIter);
RcppExport SEXP _statgenGWAS_nearestPD(SEXP xSEXP, SEXP corrSEXP, SEXP keepDiagSEXP, SEXP do2eigenSEXP, SEXP doSymSEXP, SEXP doDykstraSEXP, SEXP eigTolSEXP, SEXP convTolSEXP, SEXP posdTolSEXP, SEXP maxIterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< const bool >::type corr(corrSEXP);
    Rcpp::traits::input_parameter< const bool >::type keepDiag(keepDiagSEXP);
    Rcpp::traits::input_parameter< const bool >::type do2eigen(do2eigenSEXP);
    Rcpp::traits::input_parameter< const bool >::type doSym(doSymSEXP);
    Rcpp::traits::input_parameter< const bool >::type doDykstra(doDykstraSEXP);
    Rcpp::traits::input_parameter< const double >::type eigTol(eigTolSEXP);
    Rcpp::traits::input_parameter< const double >::type convTol(convTolSEXP);
    Rcpp::traits::input_parameter< const double >::type posdTol(posdTolSEXP);
    Rcpp::traits::input_parameter< const int >::type maxIter(maxIterSEXP);
    rcpp_result_gen = Rcpp::wrap(nearestPD(x, corr, keepDiag, do2eigen, doSym, doDykstra, eigTol, convTol, posdTol, maxIter));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_statgenGWAS_emmaEigenR", (DL_FUNC) &_statgenGWAS_emmaEigenR, 4},
    {"_statgenGWAS_emmaREMLLL", (DL_FUNC) &_statgenGWAS_emmaREMLLL, 6},
    {"_statgenGWAS_goldenSectionSearch", (DL_FUNC) &_statgenGWAS_goldenSectionSearch, 9},
    {"_statgenGWAS_emmaCPP", (DL_FUNC) &_statgenGWAS_emmaCPP, 7},
    {"_statgenGWAS_fastGLSCPP", (DL_FUNC) &_statgenGWAS_fastGLSCPP, 5},
    {"_statgenGWAS_getThr", (DL_FUNC) &_statgenGWAS_getThr, 1},
    {"_statgenGWAS_astleCPP", (DL_FUNC) &_statgenGWAS_astleCPP, 2},
    {"_statgenGWAS_IBSCPP", (DL_FUNC) &_statgenGWAS_IBSCPP, 2},
    {"_statgenGWAS_vanRadenCPP", (DL_FUNC) &_statgenGWAS_vanRadenCPP, 2},
    {"_statgenGWAS_matrixRoot", (DL_FUNC) &_statgenGWAS_matrixRoot, 1},
    {"_statgenGWAS_reduceKinship", (DL_FUNC) &_statgenGWAS_reduceKinship, 2},
    {"_statgenGWAS_nearestPD", (DL_FUNC) &_statgenGWAS_nearestPD, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_statgenGWAS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
