% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{KalmanSmoother}
\alias{KalmanSmoother}
\title{Runs a Kalman smoother}
\usage{
KalmanSmoother(F, H, R, xfT, xpT, PfT_v, PpT_v)
}
\arguments{
\item{F}{transition matrix}

\item{H}{observation matrix}

\item{R}{Observation covariance}

\item{xfT}{State estimates}

\item{PfT_v}{Variance estimates}

\item{PpT_v}{Predicted variance estimates}

\item{xpTm}{State predicted estimates}
}
\value{
List of smoothed estimates
}
\description{
Runs a Kalman smoother
}
