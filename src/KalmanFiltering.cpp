#include <RcppArmadillo.h>
#include "helper.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

//' Implementation of a Kalman filter
//' @param y Data matrix
//' @param F Observation matrix
//' @param Q State covariance
//' @param R Observation covariance
//' @param A Transition matrix
//' @param x0 Initial state vector
//' @param P0 Initial state covariance
//' @return A list with estimates
// [[Rcpp::export]]
Rcpp::List KalmanFilterCpp2(arma::mat y, arma::mat F, arma::mat Q, arma::mat R,
                            arma::mat A, arma::colvec x0, arma::mat P0) {

  const unsigned int T = y.n_rows;
  const unsigned int n = y.n_cols;
  const unsigned int r = A.n_rows;

  colvec Txittm = x0;
  mat TPttm = P0;

  double loglik = 0;
  double denom, mahal, detS = 0;
  colvec Txitt, Ires;
  mat TPtt, G, L, Icov, GG, S;

  mat xittm(T+1,r,fill::zeros);
  mat xitt(T,r,fill::zeros);
  cube Pttm(r,r,(T+1), fill::zeros);
  cube Ptt(r,r,T, fill::zeros);

  mat tF = F;
  mat tR = R;
  uvec miss;
  uvec nmiss = find_finite(F.row(0));
  uvec a(1);

  for (unsigned int t=0; t < T; t++) {
    miss = find_finite(y.row(t));
    F = tF.submat(miss, nmiss);
    R = tR.submat(miss, miss);
    a[0] = t;
    // Innovation covariance
    Icov = F * TPttm * F.t() + R;
    L = Icov.i();
    // Innovation residual
    Ires = y.submat(a, miss).t() - F * Txittm;
    // Kalman gain
    G = TPttm * F.t() * L;
    // Updated state estimate
    Txitt = Txittm + G * Ires;
    // Updated covariance estimate
    TPtt = TPttm - G * F * TPttm;
    // Compute likelihood
    S = (F * TPttm * F.t() + R).i();

    GG = F.t() * R.i() * F;
    detS = det(diagmat(R)) * det(eye(r,r) + TPttm * GG);

    denom = pow(2*datum::pi, double(n)/2.0) * sqrt(abs(detS));
    mahal = pow(norm(Ires, 2), 2);
    loglik += -0.5*mahal - log(denom);

    // Store data useful for smoothing
    xittm.row(t) = Txittm.t();
    Pttm.slice(t) = TPttm;
    xitt.row(t) = Txitt.t();
    Ptt.slice(t) = TPtt;

    // Run updates
    Txittm = A * Txitt;
    TPttm = A * TPtt * A.t() + Q;
  }

  return Rcpp::List::create(Rcpp::Named("loglik") = loglik,
                            Rcpp::Named("xitt") = xitt,
                            Rcpp::Named("xittm") = xittm,
                            Rcpp::Named("Ptt") = Ptt,
                            Rcpp::Named("Pttm") = Pttm);
}


//' Implementation of a Kalman filter
//' @param y Data matrix
//' @param H Observation matrix
//' @param Q State covariance
//' @param R Observation covariance
//' @param F Transition matrix
//' @param x0 Initial state vector
//' @param P0 Initial state covariance
// [[Rcpp::export]]
Rcpp::List KalmanFilter(arma::mat y, arma::mat H, arma::mat Q, arma::mat R,
                        arma::mat F, arma::colvec x0, arma::mat P0) {

  const unsigned int T = y.n_rows;
  const unsigned int n = y.n_cols;
  const unsigned int r = F.n_rows;

  double loglik = 0;
  mat K, Pf, Pp;
  colvec xf, xp, xe;
  // Predicted state mean and covariance
  mat xpT(T+1, r, fill::zeros);
  cube PpT(r, r, T+1, fill::zeros);

  // Filtered state mean and covariance
  mat xfT(T, r, fill::zeros);
  cube PfT(r, r, T, fill::zeros);

  mat tH = H;
  mat tR = R;
  mat S;
  uvec miss;
  uvec nmiss = find_finite(F.row(0));
  uvec a(1);

  xp = x0;
  Pp = P0;

  for (unsigned int t=0; t < T; t++) {

    // If missing observations are present at some timepoints, exclude the
    // appropriate matrix slices from the filtering procedure.
    miss = find_finite(y.row(t));
    H = tH.submat(miss, nmiss);
    R = tR.submat(miss, miss);
    a[0] = t;

    S = (H * Pp * H.t() + R).i();

    // Prediction error
    xe = y.submat(a, miss).t() - H * xp;
    // Kalman gain
    K = Pp * H.t() * S;
    // Updated state estimate
    xf = xp + K * xe;
    // Updated state covariance estimate
    Pf = Pp - K * H * Pp;

    // Compute likelihood. Skip this part if S is not positive definite.
    if (det(S) > 0) {
      loglik += -0.5 * (double(n) * log(2.0 * datum::pi) - log(det(S)) +
                        conv_to<double>::from(xe.t() * S * xe));
    }

    // Store predicted and filtered data needed for smoothing
    xpT.row(t) = xp.t();
    PpT.slice(t) = Pp;
    xfT.row(t) = xf.t();
    PfT.slice(t) = Pf;

    // Run a prediction
    xp = F * xfT.row(t).t();
    Pp = F * PfT.slice(t) * F.t() + Q;

  }

  return Rcpp::List::create(Rcpp::Named("xF") = xfT,
                            Rcpp::Named("Pf") = PfT,
                            Rcpp::Named("xP") = xpT,
                            Rcpp::Named("Pp") = PpT,
                            Rcpp::Named("loglik") = loglik);
}

//' Runs a Kalman smoother
//' @param A transition matrix
//' @param F observation matrix
//' @param R Observation covariance
//' @param xitt State estimates
//' @param xittm State predicted estimates
//' @param Ptt1 Variance estimates
//' @param Pttm1 Predicted variance estimates
//' @return List of smoothed estimates
// [[Rcpp::export]]
Rcpp::List KalmanSmootherCpp(arma::mat A, arma::mat F, arma::mat R,
                             arma::mat xitt, arma::mat xittm,
                             Rcpp::NumericVector Ptt1, Rcpp::NumericVector Pttm1) {

  const unsigned int T = xitt.n_rows;
  const unsigned int r = A.n_rows;
  const unsigned int n = F.n_rows;

  cube Ptt = array2cube(Ptt1);
  cube Pttm = array2cube(Pttm1);
  cube J(r,r,T, fill::zeros);
  cube L(n,n,T, fill::zeros);
  cube K(r,n,T, fill::zeros);

  mat xitT(T,r, fill::zeros); xitT.row(T-1) = xitt.row(T-1);
  cube PtT(r,r,T, fill::zeros); PtT.slice(T-1) = Ptt.slice(T-1);

  cube PtTm(r,r,T, fill::zeros);

  for (int t=0; t < T-1; t++) {
    J.slice(t) = Ptt.slice(t) * A.t() * Pttm.slice(t+1).i();
  }

  for (int i=0; i < T; i++) {
    L.slice(i) = (F * Pttm.slice(i) * F.t() + R).i();
    K.slice(i) = Pttm.slice(i) * F.t() * L.slice(i);
  }

  PtTm.slice(T-1) = (eye(r,r) - K.slice(T-1) * F) * A * Ptt.slice(T-2);

  for (int j=2; j < T+1; j++) {
    xitT.row(T-j) = xitt.row(T-j) + (J.slice(T-j) * (xitT.row(T+1-j)
                                                     - xittm.row(T+1-j)).t()).t();
    PtT.slice(T-j) = Ptt.slice(T-j) + J.slice(T-j) * (PtT.slice(T+1-j) -
                                      Pttm.slice(T+1-j)) * J.slice(T-j).t();
  }

  for (int j=2; j < T-1; j++) {
    PtTm.slice(T-j) = Ptt.slice(T-j) * J.slice(T-j-1).t() + J.slice(T-j)
                      * (PtTm.slice(T-j+1) - A * Ptt.slice(T-j))
                      * J.slice(T-j-1).t();
  }

  return Rcpp::List::create(Rcpp::Named("xitT")=xitT,
                            Rcpp::Named("PtT")=PtT,
                            Rcpp::Named("PtTm")=PtTm);
}

//' Runs a Kalman smoother
//' @param F transition matrix
//' @param H observation matrix
//' @param R Observation covariance
//' @param xfT State estimates
//' @param xpTm State predicted estimates
//' @param PfT_v Variance estimates
//' @param PpT_v Predicted variance estimates
//' @return List of smoothed estimates
// [[Rcpp::export]]
Rcpp::List KalmanSmoother(arma::mat F, arma::mat H, arma::mat R,
                          arma::mat xfT, arma::mat xpT,
                          Rcpp::NumericVector PfT_v, Rcpp::NumericVector PpT_v) {

  const unsigned int T = xfT.n_rows;
  const unsigned int r = F.n_rows;
  const unsigned int n = H.n_rows;

  cube PfT = array2cube(PfT_v);
  cube PpT = array2cube(PpT_v);

  cube J(r, r, T, fill::zeros);
  cube L(n, n, T, fill::zeros);
  cube K(r, n, T, fill::zeros);

  cube PsTm(r, r, T, fill::zeros);

  // Smoothed state mean and covariance
  mat xsT(T, r, fill::zeros);
  cube PsT(r,r,T, fill::zeros);
  // Initialize smoothed data with last observation of filtered data
  xsT.row(T-1) = xfT.row(T-1);
  PsT.slice(T-1) = PfT.slice(T-1);

  // cube PsTm(r,r,T, fill::zeros);
  for (int t=0; t < T-1; t++) {
    J.slice(t) = PfT.slice(t) * F.t() * PpT.slice(t+1).i();
  }

  // Smoothed state variable and covariance
  for (int j=2; j < T+1; j++) {

    xsT.row(T-j) = xfT.row(T-j) +
      (J.slice(T-j) * (xsT.row(T-j+1) - xpT.row(T-j+1)).t()).t();

    PsT.slice(T-j) = PfT.slice(T-j) +
      J.slice(T-j) * (PsT.slice(T-j+1) - PpT.slice(T-j+1)) * J.slice(T-j).t();

  }

  // Additional variables used in EM-algorithm
  for (int i=0; i < T; i++) {
    L.slice(i) = (H * PpT.slice(i) * H.t() + R).i();
    K.slice(i) = PpT.slice(i) * H.t() * L.slice(i);
  }

  PsTm.slice(T-1) = (eye(r,r) - K.slice(T-1) * H) * F * PfT.slice(T-2);

  for (int j=2; j < T-1; j++) {
    PsTm.slice(T-j) = PfT.slice(T-j) * J.slice(T-j-1).t() + J.slice(T-j)
                      * (PsTm.slice(T-j+1) - F * PfT.slice(T-j))
                      * J.slice(T-j-1).t();
  }

  return Rcpp::List::create(Rcpp::Named("xS") = xsT,
                            Rcpp::Named("Ps") = PsT,
                            Rcpp::Named("PsTm") = PsTm);
}
