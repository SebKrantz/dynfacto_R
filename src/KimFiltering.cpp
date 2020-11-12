#include <RcppArmadillo.h>
#include "helper.h"

// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;

//' Implementation of Kim filter (1994), an extension to Kalman filter
//' for dynamic linear models with Markov-switching parameters. Currently,
//' Markov switching is assumed to happen only in observation and/or
//' transition matrices.
//'
//' @param y Data matrix (\code{T x n})
//' @param R Observation equation covariance
//' @param Q State equation covariance
//' @param F1 Array of observation matrices, one matrix per state
//' @param A1 Array of transition matrices, one matrix per state
//' @param x0 Initial condition for state vector
//' @param P0 Initial condition for state covariance matrix
//' @param p Markov transition probability matrix
//' @return a list with estimates
// [[Rcpp::export]]
Rcpp::List KimFilterCpp(arma::mat y, arma::mat R, arma::mat Q,
                        Rcpp::NumericVector F1, Rcpp::NumericVector A1,
                        arma::colvec x0, arma::mat P0, arma::mat p) {

  /*
    Define all containers for further computations. Notations for variables
    and indices, where appropriate, carefully follow Kim (1994). State vector
    is denoted as 'x', its covariance as 'P'. Appended letters explicit
    whether these are updated, approximated or smoothed.

    Since computations require objects that have 4 or 5 dimensions, this is
    dealt with Field type from Armadillo library. 3-dimensional arrays are
    stored in cubes. All arrays have time as 3rd dimension (slices) while
    matrix type variables have time as 1st dimension (rows).

    Fields are 2-dimensional and indexed by states. They contain matrices
    or cubes depending on which variable is stored.
  */

  cube F = array2cube(F1);
  cube A = array2cube(A1);

  const unsigned int T = y.n_rows;
  const unsigned int n = F.n_rows;
  const unsigned int J = A.n_rows;
  const unsigned int s = p.n_rows;

  // x:   x^(i,j)_(t|t-1): predicted state vector - (2.6)
  // xU:  x^(i,j)_(t|t): updated state vector - (2.11)
  // P:   P^(i,j)_(t|t-1): predicted state covariance - (2.7)
  // Pu:  P^(i,j)_(t|t): updated state covariance - (2.12)
  // eta: eta^(i,j)_(t|t-1): conditional forecast error - (2.8)
  // H:   H^(i,j)_(t): conditional variance of forecast error - (2.9)
  // K:   K^(i,j)_(t): Kalman gain - (2.10)
  field<mat> x(s, s);       // T x J
  field<mat> xU(s, s);      // T x J
  field<cube> P(s, s);      // J x J x T
  field<cube> Pu(s, s);     // J x J x T
  field<mat> eta(s, s);     // T x n
  field<cube> H(s, s);      // n x n x T
  field<cube> K(s, s);      // J x n x T

  // lik: f(y_t, S_{t-1} = i, S_t = j | t-1): joint conditional density - (2.16)
  // loglik: log of (2.16)
  cube lik(s, s, T, fill::zeros);
  cube loglik(s, s, T, fill::zeros);

  // Pr[S_(t-1) = i, S_t = j | t-1 ]: (2.15)
  // Pr[S_(t-1) = i, S_t = j | t ]: (2.17)
  // Pr[S_t = j | t-1 ]: used only for the smoothing part
  // Pr[S_t = j | t ]: (2.18)
  cube jointP_fut(s, s, T, fill::zeros);
  cube jointP_cur(s, s, T, fill::zeros);
  mat stateP_fut(T, s, fill::zeros);
  mat stateP(T, s, fill::zeros);

  // x^(j)_(t|t): approximate state vector conditional on S_j - (2.13)
  // P^(j)_(t|t): approximate state covariance conditional on S_j - (2.14)
  cube xA(J, s, T, fill::zeros);     // J x s x T
  field<cube> Pa(s);                 // J x J x T
  colvec result(T, fill::zeros);     // T x 1
  // End of definitions. Initialize all fields with fixed size matrices.
  mat mTJ(T, J, fill::zeros), mTn(T,n, fill::zeros);
  cube cJJT(J, J, T, fill::zeros), cnnT(n,n,T, fill::zeros), cJnT(J,n,T, fill::zeros);
  for (unsigned int i=0; i < s; i++) {
    for (unsigned int j=0; j < s; j++) {
      x(i,j) = mTJ; xU(i,j) = mTJ; eta(i,j) = mTn;
      P(i,j) = cJJT; Pu(i,j) = cJJT; H(i,j) = cnnT; K(i,j) = cJnT;
    }
    Pa(i) = cJJT;
  }

  double detH = 0, etaH = 0;

  // Some initial conditions to get started
  jointP_fut.slice(0) = "0.25 0.25; 0.25 0.25";
  for (unsigned int i=0; i < s; i++) {
    for (unsigned int j=0; j < s; j++) {
      x(i,j).row(0) = x0.t();
      P(i,j).slice(0) = P0;
    }
  }

  for (unsigned int t=0; t < T; t++) {
    for (unsigned int j=0; j < s; j++) {
      for (unsigned int i=0; i < s; i++) {

        if (t > 0) {
          x(i,j).row(t) = (A.slice(j) * xA.slice(t-1).col(i)).t();
          P(i,j).slice(t) = A.slice(j) * Pa(i).slice(t-1) * A.slice(j).t() + Q;
          jointP_fut(i,j,t) = p(i,j) * accu(jointP_cur.slice(t-1).col(i));
        }

        eta(i,j).row(t) = y.row(t) - (F.slice(j) * x(i,j).row(t).t()).t();
        H(i,j).slice(t) = F.slice(j) * P(i,j).slice(t) * F.slice(j).t() + R;
        K(i,j).slice(t) = P(i,j).slice(t) * F.slice(j).t() * H(i,j).slice(t).i();

        xU(i,j).row(t) = x(i,j).row(t) + (K(i,j).slice(t) * eta(i,j).row(t).t()).t();
        Pu(i,j).slice(t) = (eye(J,J) - K(i,j).slice(t) * F.slice(j)) * P(i,j).slice(t);
        // Computing (log)-likelihood
        detH = det(H(i,j).slice(t));
        etaH = as_scalar(eta(i,j).row(t) * H(i,j).slice(t).i() * eta(i,j).row(t).t());
        loglik(i,j,t) = -0.5 * (double(n) * log(2.0 * datum::pi) + log(detH) + etaH) +
          log(jointP_fut(i,j,t));
        lik(i,j,t) = exp(loglik(i,j,t));

        jointP_cur(i,j,t) = lik(i,j,t);
      }
      // Technically, there should be accu(lik.slice(t)) term but it cancels out
      // and is computed later
      stateP(t,j) = accu(jointP_cur.slice(t).col(j));
      stateP_fut(t,j) = accu(jointP_fut.slice(t).col(j));
      // Compute probability-filtered state process and its covariance
      for (int i=0; i < s; i++) {
        xA.slice(t).col(j) += (xU(i,j).row(t) * jointP_cur(i,j,t) / stateP(t,j)).t();
      }
      for (int i=0; i < s; i++) {
        Pa(j).slice(t) += (Pu(i,j).slice(t)
                           + (xA.slice(t).col(j) - xU(i,j).row(t).t())
                           * (xA.slice(t).col(j) - xU(i,j).row(t).t()).t())
          * (jointP_cur(i,j,t) / stateP(t,j));
      }
    }
    jointP_cur.slice(t) = jointP_cur.slice(t) / accu(lik.slice(t));
    stateP.row(t) = stateP.row(t) / accu(lik.slice(t));
    result(t) = log(accu(lik.slice(t)));

  }

  return Rcpp::List::create(Rcpp::Named("loglik") = loglik,
                            Rcpp::Named("result") = result,
                            Rcpp::Named("xA") = xA,
                            Rcpp::Named("lik") = loglik,
                            Rcpp::Named("Pa") = Pa,
                            Rcpp::Named("x") = x,
                            Rcpp::Named("P") = P,
                            Rcpp::Named("stateP") = stateP,
                            Rcpp::Named("stateP_fut") = stateP_fut);
}

// [[Rcpp::export]]
Rcpp::List KimSmoother(Rcpp::NumericVector F1,
                       Rcpp::NumericVector xA1,
                       Rcpp::NumericVector Pa1,
                       Rcpp::NumericVector x1,
                       Rcpp::NumericVector P1,
                       arma::mat stateP, arma::mat stateP_fut, arma::mat p) {

  cube F = array2cube(F1);
  cube xA = array2cube(xA1);
  field<cube> Pa = array2field1cube(Pa1);
  field<mat> x = array2field2mat(x1);
  field<cube> P = array2field2cube(P1);

  // Define all containers for further computations. Notations for variables and
  // indices, where appropriate, carefully follow Kim (1994). State vector is
  // denoted as 'x', its covariance as 'P'. Appended letters explicit whether
  // these are updated, approximated or smoothed.
  const unsigned int T = stateP.n_rows;
  const unsigned int J = xA.n_cols;
  const unsigned int s = p.n_rows;

  // Pr[S_t = j, S_(t+1) = k | T]: (2.20)
  // Pr[S_t = j | T]: (2.21)
  cube jointPs(s, s, T, fill::zeros);
  mat ProbS(T, s, fill::zeros);

  // xS: x^(j,k)_(t|T): inference of x_t based on full sample - (2.24)
  // Ps: P^(j,k)_(t|T): covariance matrix of x^(j,k)_(t|T) - (2.25)
  // Ptilde: helper matrix as defined after (2.25)
  field<mat> xS(s, s);         // T x J
  field<cube> Ps(s, s);        // J x J x T
  field<cube> Ptilde(s, s);     // J x J x T

  // xAS: x^(j)_(t|T): smoothed and approximated state vector conditional on S_j (2.26)
  // Pas: P^(j)_(t|T): smoothed and approximated state covariance conditional on S_j (2.27)
  // xF: x_(t|T): state-weighted [F]inal state vector (2.28)
  // Pf: P_(t|T): state-weighted [f]inal state covariance
  cube xAS(J, s, T, fill::zeros);     // J x s x T
  field<cube> Pas(s);                 // J x J x T
  mat xF(T, J, fill::zeros);
  cube Pf(s, s, T, fill::zeros);

  for (unsigned int i=0; i < s; i++) {
    for (unsigned int j=0; j < s; j++) {
      xS(i,j) = mat(T, J, fill::zeros);
      Ps(i,j) = cube(J, J, T, fill::zeros);
      Ptilde(i,j) = cube(J, J, T, fill::zeros);
    }
    Pas(i) = cube(J, J, T, fill::zeros);
    // Initial condition for smoothing covariance
    Pas(i).slice(T-1) = Pa(i).slice(T-1);
  }

  // Initial conditions for smoothing variables
  ProbS.row(T-1) = stateP.row(T-1);
  xAS.slice(T-1) = xA.slice(T-1);

  // Construct probability objects
  for (int t=(T-2); t >= 0; t--) {
    for (unsigned int j=0; j < s; j++) {
      for (unsigned int k=0; k < s; k++) {
        jointPs(j, k, t) = ProbS(t+1, k) * stateP(t,j) * p(j,k) / stateP_fut(t+1, k);
      }
      ProbS(t, j) = accu(jointPs.slice(t).row(j));
    }
  }

  for (int t=(T-2); t >= 0; t--) {
    for (unsigned int j=0; j < s; j++) {
      for (unsigned int k=0; k < s; k++) {

        Ptilde(j, k).slice(t) = Pa(j).slice(t) * F.slice(k).t() * P(j,k).slice(t+1).i();

        xS(j,k).row(t) = xA.slice(t).col(j).t() +
          (Ptilde(j,k).slice(t) * (xAS.slice(t+1).col(k) - x(j,k).row(t+1).t())).t();
        Ps(j,k).slice(t) = Pa(j).slice(t) + Ptilde(j,k).slice(t) *
          (Pas(k).slice(t+1) - P(j,k).slice(t+1)) * Ptilde(j,k).slice(t).t();

        xAS.slice(t).col(j) = xAS.slice(t).col(j) + jointPs(j,k,t) * xS(j,k).row(t).t();
      }

      xAS.slice(t).col(j) = xAS.slice(t).col(j) / ProbS(t,j);
      // xAS needs to be calculated first before using it for Pas
      for (int k = 0; k < s; k++) {
        Pas(j).slice(t) = Pas(j).slice(t) + jointPs(j,k,t) *
          (Ps(j, k).slice(t) + (xAS.slice(t).col(j) - xS(j,k).row(t).t()) *
           (xAS.slice(t).col(j).t() - xS(j,k).row(t)));
      }
      Pas(j).slice(t) = Pas(j).slice(t) / ProbS(t,j);

    }
  }

  for (unsigned int t = 0; t < T; t++) {
    for (unsigned int j = 0; j < J; j++) {
      xF.row(t) = xF.row(t) + xAS.slice(t).col(j).t() * ProbS(t,j);
      Pf.slice(t) = Pf.slice(t) + Pas(j).slice(t) * ProbS(t,j);
    }
  }

  return Rcpp::List::create(Rcpp::Named("xF") = xF,
                            Rcpp::Named("Pf") = Pf,
                            Rcpp::Named("ProbS") = ProbS);

}
