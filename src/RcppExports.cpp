// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// EstepCpp
Rcpp::List EstepCpp(arma::mat y, arma::mat H, arma::mat Q, arma::mat R, arma::mat F, arma::colvec x0, arma::mat P0);
RcppExport SEXP dynfactoR_EstepCpp(SEXP ySEXP, SEXP HSEXP, SEXP QSEXP, SEXP RSEXP, SEXP FSEXP, SEXP x0SEXP, SEXP P0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type H(HSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q(QSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type R(RSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type F(FSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P0(P0SEXP);
    rcpp_result_gen = Rcpp::wrap(EstepCpp(y, H, Q, R, F, x0, P0));
    return rcpp_result_gen;
END_RCPP
}
// KalmanFilterCpp2
Rcpp::List KalmanFilterCpp2(arma::mat y, arma::mat F, arma::mat Q, arma::mat R, arma::mat A, arma::colvec x0, arma::mat P0);
RcppExport SEXP dynfactoR_KalmanFilterCpp2(SEXP ySEXP, SEXP FSEXP, SEXP QSEXP, SEXP RSEXP, SEXP ASEXP, SEXP x0SEXP, SEXP P0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type F(FSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q(QSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type R(RSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P0(P0SEXP);
    rcpp_result_gen = Rcpp::wrap(KalmanFilterCpp2(y, F, Q, R, A, x0, P0));
    return rcpp_result_gen;
END_RCPP
}
// KalmanFilter
Rcpp::List KalmanFilter(arma::mat y, arma::mat H, arma::mat Q, arma::mat R, arma::mat F, arma::colvec x0, arma::mat P0);
RcppExport SEXP dynfactoR_KalmanFilter(SEXP ySEXP, SEXP HSEXP, SEXP QSEXP, SEXP RSEXP, SEXP FSEXP, SEXP x0SEXP, SEXP P0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type H(HSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q(QSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type R(RSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type F(FSEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P0(P0SEXP);
    rcpp_result_gen = Rcpp::wrap(KalmanFilter(y, H, Q, R, F, x0, P0));
    return rcpp_result_gen;
END_RCPP
}
// KalmanSmootherCpp
Rcpp::List KalmanSmootherCpp(arma::mat A, arma::mat F, arma::mat R, arma::mat xitt, arma::mat xittm, Rcpp::NumericVector Ptt1, Rcpp::NumericVector Pttm1);
RcppExport SEXP dynfactoR_KalmanSmootherCpp(SEXP ASEXP, SEXP FSEXP, SEXP RSEXP, SEXP xittSEXP, SEXP xittmSEXP, SEXP Ptt1SEXP, SEXP Pttm1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type F(FSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type R(RSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xitt(xittSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xittm(xittmSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Ptt1(Ptt1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Pttm1(Pttm1SEXP);
    rcpp_result_gen = Rcpp::wrap(KalmanSmootherCpp(A, F, R, xitt, xittm, Ptt1, Pttm1));
    return rcpp_result_gen;
END_RCPP
}
// KalmanSmoother
Rcpp::List KalmanSmoother(arma::mat F, arma::mat H, arma::mat R, arma::mat xfT, arma::mat xpT, Rcpp::NumericVector PfT_v, Rcpp::NumericVector PpT_v);
RcppExport SEXP dynfactoR_KalmanSmoother(SEXP FSEXP, SEXP HSEXP, SEXP RSEXP, SEXP xfTSEXP, SEXP xpTSEXP, SEXP PfT_vSEXP, SEXP PpT_vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type F(FSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type H(HSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type R(RSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xfT(xfTSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xpT(xpTSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type PfT_v(PfT_vSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type PpT_v(PpT_vSEXP);
    rcpp_result_gen = Rcpp::wrap(KalmanSmoother(F, H, R, xfT, xpT, PfT_v, PpT_v));
    return rcpp_result_gen;
END_RCPP
}
// KimFilterCpp
Rcpp::List KimFilterCpp(arma::mat y, arma::mat R, arma::mat Q, Rcpp::NumericVector F1, Rcpp::NumericVector A1, arma::colvec x0, arma::mat P0, arma::mat p);
RcppExport SEXP dynfactoR_KimFilterCpp(SEXP ySEXP, SEXP RSEXP, SEXP QSEXP, SEXP F1SEXP, SEXP A1SEXP, SEXP x0SEXP, SEXP P0SEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type R(RSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Q(QSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type F1(F1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type A1(A1SEXP);
    Rcpp::traits::input_parameter< arma::colvec >::type x0(x0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type P0(P0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(KimFilterCpp(y, R, Q, F1, A1, x0, P0, p));
    return rcpp_result_gen;
END_RCPP
}
// KimSmoother
Rcpp::List KimSmoother(Rcpp::NumericVector F1, Rcpp::NumericVector xA1, Rcpp::NumericVector Pa1, Rcpp::NumericVector x1, Rcpp::NumericVector P1, arma::mat stateP, arma::mat stateP_fut, arma::mat p);
RcppExport SEXP dynfactoR_KimSmoother(SEXP F1SEXP, SEXP xA1SEXP, SEXP Pa1SEXP, SEXP x1SEXP, SEXP P1SEXP, SEXP statePSEXP, SEXP stateP_futSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type F1(F1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type xA1(xA1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Pa1(Pa1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x1(x1SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type P1(P1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type stateP(statePSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type stateP_fut(stateP_futSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(KimSmoother(F1, xA1, Pa1, x1, P1, stateP, stateP_fut, p));
    return rcpp_result_gen;
END_RCPP
}