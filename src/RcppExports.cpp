// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// calc_TMTI_inf
List calc_TMTI_inf(NumericVector pvals);
RcppExport SEXP _TMTI_calc_TMTI_inf(SEXP pvalsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type pvals(pvalsSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_TMTI_inf(pvals));
    return rcpp_result_gen;
END_RCPP
}
// calc_truncatedTMTI_inf
List calc_truncatedTMTI_inf(NumericVector pvals, int m_full);
RcppExport SEXP _TMTI_calc_truncatedTMTI_inf(SEXP pvalsSEXP, SEXP m_fullSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type pvals(pvalsSEXP);
    Rcpp::traits::input_parameter< int >::type m_full(m_fullSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_truncatedTMTI_inf(pvals, m_full));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_TMTI_calc_TMTI_inf", (DL_FUNC) &_TMTI_calc_TMTI_inf, 1},
    {"_TMTI_calc_truncatedTMTI_inf", (DL_FUNC) &_TMTI_calc_truncatedTMTI_inf, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_TMTI(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
