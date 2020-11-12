#' Estimates a dynamic factor model based on Doz, Gianone & Reichlin (2011)
#'
#' @param data Data matrix with time in rows and variables in columns.
#' @param r Number of static factors to estimate.
#' @param p Lag order for factors. It is assumed that factors follow VAR(p)
#'   model.
#' @param q Number of dynamic factors, must be equal to or less than r. Dynamic
#'   factors refer essentially to the number of principal components relevant
#'   for the state covariance estimation.
#' @param max_iter Maximum number of iterations in the EM-algorithm.
#' @param threshold Threshold for algorithm convergence
#' @param lower_set In order to keep an invertible observation covariance
#'   matrix, diagonal values are constrained to be such that
#'   \code{diag(R)[diag(R) < lower_set] <- lower_set}.
#' @param rQ Restrictions on system state covariance. Currently, either no
#'   restrictions can be set or \code{Q} can be fixed as identity.
#' @param rC Restrictions on factor loading matrix. Currently, either no
#'   restrictions can be set or upper matrix of \code{C} can be set to 0.
#' @param ... Further arguments that are currently unused
#' @return 3 types of factor estimates, namely principal component estimate, two
#'   step estimate based on PCA and Kalman filtering and QML estimate based on
#'   EM-algorithm. PCA estimator is not able to deal with missing data.
#' @importFrom MASS ginv
#' @importFrom stats cov na.omit sd
#' @export
#' @examples
#' \dontrun{
#' set.seed(314)
#' x <- matrix(rnorm(50*10), 50, 10)
#' W <- as.logical(matrix(rbinom(50*10, 1, 0.1), 50, 10))
#' x[W] <- NA
#' dfm(x, 2, 2, 1)
#' }
dfm <- function(data, r, p, q = r, max_iter = 100,
                threshold = 1e-4, lower_set = 1e-7,
                rQ = c("none", "identity"),
                rC = c("none", "upper"), ...) {

  rQ <- match.arg(rQ)
  rC <- match.arg(rC)

  ## C - Observation matrix
  ## A - State transition matrix
  ## Q - State covariance matrix
  ## R - Observation covariance matrix

  if (q > r) stop("r must be larger than q.")

  T <- dim(data)[1]
  n <- dim(data)[2]

  ## Data is stantardized so that PCA can work properly. Mean and variance of
  ## the original data are saved added back at the end.
  x <- apply(data, 2, function(z) {
    (z - mean(z, na.rm = TRUE)) / sd(z, na.rm = TRUE)
  })
  Mx <- apply(data, 2, mean, na.rm = TRUE)
  Wx <- apply(data, 2, sd, na.rm = TRUE)
  W <- !is.na(x)

  ## State transition matrix A consists of two parts. In particular, the upper
  ## dynamic part and the lower invariable identity part to ensure time lag
  ## coherence.
  A <- rbind(matrix(0, nrow = r, ncol = r*p),
             diag(1, nrow = r*(p-1), ncol = r*p))

  ## State covariance matrix Q
  Q <- matrix(0, nrow = p*r, ncol = p*r)
  Q[1:r, 1:r] <- diag(1, r)

  eigen.decomp <- eigen(cov(x, use="complete.obs"))
  v <- eigen.decomp$vectors[,1:r]
  d <- eigen.decomp$values[1:r]

  chi <- x %*% v %*% t(v)
  d <- diag(1, r)
  F <- x %*% v

  ## PCA solution refers to first r principal components
  F_pc <- F

  F <- na.omit(F)

  ## If there are any dynamic factors, VAR(p) model is estimated
  ## to initialize their parameters.
  if (p > 0) {
    ## ML estimator for VAR(p) model when Q is restricted
    fit <- VAR(F, p)
    if (rQ == "identity") {
      A[1:r, 1:(r*p)] <- t(fit$A)
      Q[1:r, 1:r] <- diag(1, r)
    } else {
      A[1:r, 1:(r*p)] <- t(fit$A)
      H <- cov(fit$res)

      ## This only extracts the variance explained by the dynamic components
      if (r > q) {
        q.decomp <- eigen(H)
        P <- q.decomp$vectors[,1:q, drop=FALSE]
        M <- q.decomp$values[1:q]
        if (q == 1) {
          P <- P * P[1,]
          Q[1:r, 1:r] <- P %*% t(P) * M
        } else {
          P <- P %*% diag(sign(P[1,]))
          Q[1:r, 1:r] <- P %*% diag(M) %*% t(P)
        }
      } else {
        Q[1:r, 1:r] <- H
      }
    }
  }

  ## Observation covariance matrix R
  R <- diag(diag(cov(x - chi, use = "complete.obs")))

  ## Initial values for state mean and covariance
  x0 <- fit$X[1,]
  P0 <- matrix(ginv(kronecker(A,A)) %*% as.numeric(Q),
               ncol = r*p, nrow = r*p)

  ## Observation matrix C initialized with eigenvectors
  C <- cbind(v, matrix(0, nrow = n, ncol = r*(p-1)))

  ## Run standartized data through Kalman filter and smoother once
  kf_res <- KalmanFilter(x, C, Q, R, A, x0, P0)
  ks_res <- with(kf_res, KalmanSmoother(A, C, R, xF, xP, Pf, Pp))

  ## Two-step solution is state mean from the Kalman smoother
  F_kal <- ks_res$xS


  if (rC == "upper" & (r > 1)) {
    dimC <- dim(C[, 1:r])
    rK <- rep(0, (r-1)*r/2)
    irC <- which(matrix(upper.tri(C[, 1:r]) + 0) == 1)
    rH <- matrix(0, nrow = length(rK), ncol = prod(dimC))
    for (i in 1:length(rK)) {
      rH[i, irC[i]] <- 1
    }
  }

  previous_loglik <- -.Machine$double.xmax
  num_iter <- 0
  converged <- FALSE

  while ((num_iter < max_iter) & !converged) {

    ## E-step will return a list of sufficient statistics, namely second
    ## (cross)-moments for latent and observed data. This is then plugged back
    ## into M-step.
    em_res <- EstepCpp(x, C, Q, R, A, x0, P0)
    beta <- em_res$beta_t
    gamma <- em_res$gamma_t
    delta <- em_res$delta_t
    gamma1 <- em_res$gamma1_t
    gamma2 <- em_res$gamma2_t
    P1sum <- em_res$V1 + em_res$x1 %*% t(em_res$x1)
    x1sum <- em_res$x1
    loglik <- em_res$loglik_t

    num_iter <- num_iter + 1

    ## M-step computes model parameters as a function of the sufficient
    ## statistics that were computed with the E-step. Iterate the procedure
    ## until convergence. Due to the model specification, likelihood maximiation
    ## in the M-step is just an OLS estimation. In particular, X_t = C*F_t and
    ## F_t = A*F_(t-1).

    if (rC == "upper" & (r > 1)) {

      fp <- matrix(delta[, 1:r] %*% ginv(gamma[1:r, 1:r]))
      kronCR <- kronecker(ginv(gamma[1:r,1:r]), R)
      sp <- kronCR %*% t(rH) %*% ginv(rH %*% kronCR %*% t(rH)) %*% (rK - rH %*% fp)
      C[, 1:r] <- matrix(fp + sp, nrow = dimC[1], ncol = dimC[2])

    } else {
      C[, 1:r] <- delta[, 1:r] %*% ginv(gamma[1:r, 1:r])
    }

    if (p > 0) {

      A_update <- beta[1:r, 1:(r*p), drop = FALSE] %*% solve(gamma1[1:(r*p), 1:(r*p)])
      A[1:r, 1:(r*p)] <- A_update

      if (rQ != "identity") {
        H <- (gamma2[1:r, 1:r] - A_update %*% t(beta[1:r, 1:(r*p), drop=FALSE])) / (T-1)
        if (r > q) {
          h.decomp <- svd(H)
          P <- h.decomp$v[, 1:q, drop=FALSE]
          M <- h.decomp$d[1:q]
          if (q == 1) {
            P <- P * P[1,]
            Q[1:r, 1:r] <- P %*% t(P) * M
          } else {
            P <- P %*% diag(sign(P[1,]))
            Q[1:r, 1:r] <- P %*% diag(M) %*% t(P)
          }
        } else {
          Q[1:r, 1:r] <- H
        }
      }
    }

    xx <- as.matrix(na.omit(x))
    R <- (t(xx) %*% xx - C %*% t(delta)) / T
    RR <- diag(R); RR[RR < lower_set] <- lower_set; R <- diag(RR)

    R <- diag(diag(R))

    ## Assign new initial values for next EM-algorithm step
    x0 <- x1sum
    P0 <- P1sum - x0 %*% t(x0)

    converged <- em_converged(loglik, previous_loglik,
                              threshold = threshold)
    previous_loglik <- loglik

    ## Iterate at least 25 times
    if (num_iter < 25) converged <- FALSE

  }

  if (converged == TRUE)
    message("Converged after ", num_iter, " iterations.")
  else
    warning("Maximum number of iterations reached.")

  ## Run the Kalman filtering and smoothing step for the last time
  ## with optimal estimates
  kf <- KalmanFilter(x, C, Q, R, A, x0, P0)
  ks <- KalmanSmoother(A, C, R, kf$xF, kf$xP, kf$Pf, kf$Pp)

  xS <- ks$xS

#  chi <- t(xsmooth) %*% t(C) %*% diag(Wx) + kronecker(matrix(1,T,1), t(Mx))
#  F_hat <- t(xsmooth[1:r,, drop=FALSE])
  F_hat <- xS
  final_object <- list("pca" = F_pc,
                       "twostep" = F_kal,
                       "qml" = F_hat,
                       "A" = A[1:r, ],
                       "C" = C[, 1:r],
                       "Q" = Q[1:q, 1:q],
                       "R" = R)
  class(final_object) <- append("dfm", class(final_object))

  final_object

}
