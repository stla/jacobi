// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "jacobi_types.h"
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// modulo
double modulo(double a, double p);
RcppExport SEXP _jacobi_modulo(SEXP aSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(modulo(a, p));
    return rcpp_result_gen;
END_RCPP
}
// ljtheta2_cpp
cplx ljtheta2_cpp(cplx z, cplx tau);
RcppExport SEXP _jacobi_ljtheta2_cpp(SEXP zSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< cplx >::type z(zSEXP);
    Rcpp::traits::input_parameter< cplx >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(ljtheta2_cpp(z, tau));
    return rcpp_result_gen;
END_RCPP
}
// jtheta2_cpp
cplx jtheta2_cpp(cplx z, cplx tau);
RcppExport SEXP _jacobi_jtheta2_cpp(SEXP zSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< cplx >::type z(zSEXP);
    Rcpp::traits::input_parameter< cplx >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(jtheta2_cpp(z, tau));
    return rcpp_result_gen;
END_RCPP
}
// ljtheta1_cpp
cplx ljtheta1_cpp(cplx z, cplx tau);
RcppExport SEXP _jacobi_ljtheta1_cpp(SEXP zSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< cplx >::type z(zSEXP);
    Rcpp::traits::input_parameter< cplx >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(ljtheta1_cpp(z, tau));
    return rcpp_result_gen;
END_RCPP
}
// jtheta1_cpp
cplx jtheta1_cpp(cplx z, cplx tau);
RcppExport SEXP _jacobi_jtheta1_cpp(SEXP zSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< cplx >::type z(zSEXP);
    Rcpp::traits::input_parameter< cplx >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(jtheta1_cpp(z, tau));
    return rcpp_result_gen;
END_RCPP
}
// ljtheta3_cpp
cplx ljtheta3_cpp(cplx z, cplx tau);
RcppExport SEXP _jacobi_ljtheta3_cpp(SEXP zSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< cplx >::type z(zSEXP);
    Rcpp::traits::input_parameter< cplx >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(ljtheta3_cpp(z, tau));
    return rcpp_result_gen;
END_RCPP
}
// jtheta3_cpp
cplx jtheta3_cpp(cplx z, cplx tau);
RcppExport SEXP _jacobi_jtheta3_cpp(SEXP zSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< cplx >::type z(zSEXP);
    Rcpp::traits::input_parameter< cplx >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(jtheta3_cpp(z, tau));
    return rcpp_result_gen;
END_RCPP
}
// ljtheta4_cpp
cplx ljtheta4_cpp(cplx z, cplx tau);
RcppExport SEXP _jacobi_ljtheta4_cpp(SEXP zSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< cplx >::type z(zSEXP);
    Rcpp::traits::input_parameter< cplx >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(ljtheta4_cpp(z, tau));
    return rcpp_result_gen;
END_RCPP
}
// jtheta4_cpp
cplx jtheta4_cpp(cplx z, cplx tau);
RcppExport SEXP _jacobi_jtheta4_cpp(SEXP zSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< cplx >::type z(zSEXP);
    Rcpp::traits::input_parameter< cplx >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(jtheta4_cpp(z, tau));
    return rcpp_result_gen;
END_RCPP
}
// Image_eta
Rcpp::CharacterMatrix Image_eta(Rcpp::NumericVector x, cplx gamma, double t);
RcppExport SEXP _jacobi_Image_eta(SEXP xSEXP, SEXP gammaSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< cplx >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(Image_eta(x, gamma, t));
    return rcpp_result_gen;
END_RCPP
}
// Image_E4
Rcpp::CharacterMatrix Image_E4(Rcpp::NumericVector x, cplx gamma, double t);
RcppExport SEXP _jacobi_Image_E4(SEXP xSEXP, SEXP gammaSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< cplx >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(Image_E4(x, gamma, t));
    return rcpp_result_gen;
END_RCPP
}
// Image_E6
Rcpp::CharacterMatrix Image_E6(Rcpp::NumericVector x, cplx gamma, double t);
RcppExport SEXP _jacobi_Image_E6(SEXP xSEXP, SEXP gammaSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< cplx >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(Image_E6(x, gamma, t));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_jacobi_modulo", (DL_FUNC) &_jacobi_modulo, 2},
    {"_jacobi_ljtheta2_cpp", (DL_FUNC) &_jacobi_ljtheta2_cpp, 2},
    {"_jacobi_jtheta2_cpp", (DL_FUNC) &_jacobi_jtheta2_cpp, 2},
    {"_jacobi_ljtheta1_cpp", (DL_FUNC) &_jacobi_ljtheta1_cpp, 2},
    {"_jacobi_jtheta1_cpp", (DL_FUNC) &_jacobi_jtheta1_cpp, 2},
    {"_jacobi_ljtheta3_cpp", (DL_FUNC) &_jacobi_ljtheta3_cpp, 2},
    {"_jacobi_jtheta3_cpp", (DL_FUNC) &_jacobi_jtheta3_cpp, 2},
    {"_jacobi_ljtheta4_cpp", (DL_FUNC) &_jacobi_ljtheta4_cpp, 2},
    {"_jacobi_jtheta4_cpp", (DL_FUNC) &_jacobi_jtheta4_cpp, 2},
    {"_jacobi_Image_eta", (DL_FUNC) &_jacobi_Image_eta, 3},
    {"_jacobi_Image_E4", (DL_FUNC) &_jacobi_Image_E4, 3},
    {"_jacobi_Image_E6", (DL_FUNC) &_jacobi_Image_E6, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_jacobi(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
