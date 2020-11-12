#' Dynamic factor model with Markov-switching states
#'
#' This function extends in a sense \code{\link{dfm}} in that it allows
#' observation or transition matrices to follow a Markov-switching
#' process, so that they are state-dependent.
#'
#' Note that this method does not implement Chang-Jin Kim algorithm. It is based
#' on the same idea as the \code{\link{dfm}} function with the adaptation that
#' Kalman filter is replaced by Kim filter which is able to deal with multiple
#' states. Optimization is done over a list of parameters by a simple call to
#' \code{\link{optim}} instead of an EM-algorithm.
#'
#' Currently, state equation covariance matrix is restricted to identity.
#'
#' @param data Data matrix
#' @param nf Number of factors
#' @param ns Number of states
#' @param x0 Initial value for state vector
#' @param P0 Initial value for state covariance, i.e. uncertainty of the initial
#'   state value vector
#' @param init List with initial values for A, F, R, p. If not supplied,
#'   \code{\link{dfm}} function will be run and the estimates will be used as
#'   initial values.
#' @importFrom stats optim runif sd
#' @importFrom utils head tail
#' @export
dfmMS <- function(data, nf = 1, ns = 2,
                  x0 = rep(0, nf), P0 = diag(1, nf),
                  init, ...) {

  if (any(is.na(data)))
    stop("Missing observations are currently not allowed.")

  ## n = number of observations
  n <- ncol(data)

  ## Normalize data
  x <- apply(data, 2, function(z) {
    (z - mean(z, na.rm = TRUE)) / sd(z, na.rm = TRUE)
  })

  if (missing(init)) {
    init <- list()
    init$A <- array(0, c(nf, nf, ns))
    init$F <- array(1, c(n, nf, ns))
    init$R <- rep(1, n)
    init$p <- c(0.5, 0.5)
  }

  ## State equation covariance is restricted to identity
  Q <- diag(1, nf)

  ## Setup the dimensions
  dimA <- nf * nf * ns
  dimF <- n * nf * ns
  dimAF <- dimA + dimF
  dimAFR <- dimAF + n

  KimFilterOptim <- function(pars) {

    A <- array(head(pars, dimA), c(nf, nf, ns))
    F <- array(pars[(dimA+1):(dimAF)], c(n, nf, ns))
    R <- diag(pars[(dimAF+1):(dimAFR)])

    dp <- tail(pars, ns); r <- 1 - dp
    p <- diag(dp) + r - diag(r)

    ## Return log-likelihood estimate
    sum(KimFilter(x, R, Q, F, A, x0, P0, p)$loglik)
  }

  optimP <- optim(par = unlist(init),
                  fn = KimFilterOptim, ...)

  ## Reorganize parameters
  pars <- optimP$par
  A_hat <- array(pars[1:dimA], c(nf, nf, ns))
  F_hat <- array(pars[(dimA+1):dimAF], c(n, nf, ns))
  R_hat <- diag(pars[(dimAF+1):dimAFR])
  p_hat <- tail(pars, ns); p_hat <- diag(p_hat) + (1-p_hat) - diag(1-p_hat)

  ## Run a final filtering and smoothing step
  kf <- KimFilter(x, R_hat, Q, F_hat, A_hat, x0, P0, p_hat)
  ks <- KimSmoother(A_hat, kf$xA, kf$Pa, kf$x, kf$P, kf$stateP, kf$stateP_fut, p_hat)

  list("A" = A_hat, "F" = F_hat, "R" = R_hat,
       "p" = p_hat, "xF" = ks$xF,
       "Pf" = ks$Pf, "ProbS" = ks$ProbS)
}
