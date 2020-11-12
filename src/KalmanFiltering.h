Rcpp::List KalmanFilter(arma::mat y, arma::mat F, arma::mat Q, arma::mat R,
                        arma::mat A, arma::colvec x0, arma::mat P0);

Rcpp::List KalmanSmoother(arma::mat A, arma::mat F, arma::mat R,
                          arma::mat xitt, arma::mat xittm,
                          Rcpp::NumericVector Ptt1, Rcpp::NumericVector Pttm1);

Rcpp::List KalmanFilterCpp2(arma::mat y, arma::mat F, arma::mat Q, arma::mat R,
                            arma::mat A, arma::colvec x0, arma::mat P0);

Rcpp::List KalmanSmootherCpp(arma::mat A, arma::mat F, arma::mat R,
                             arma::mat xitt, arma::mat xittm,
                             Rcpp::NumericVector Ptt1, Rcpp::NumericVector Pttm1);

