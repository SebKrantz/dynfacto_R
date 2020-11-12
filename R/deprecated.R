#' Computation of the expectation step in the EM-algorithm.

#' @param y Data matrix
#' @param A System state matrix
#' @param C Observation matrix
#' @param R Observation equation variance
#' @param Q System state equation variance
#' @param W Logical matrix with dim(W) = dim(y)
#' indicating missing observations
#' @param initx Initial value for state variable
#' @param initV Initial value for state matrix
#' @return Sufficient statistics used for M-step
Estep <- function(y, A, C, Q, R, initx, initV, W) {

    os <- dim(y)[1]
    T <- dim(y)[2]
    ss <- nrow(A)

    kf <- K_filter(initx, initV, t(y), A, C, R, Q)
    ks <- K_smoother(A, kf$xitt, kf$xittm, kf$Ptt, kf$Pttm, C, R, W)

  #  kf <- KalmanFilterCpp(t(y), C, Q, R, A, initx, initV)
  #  ks <- KalmanSmootherCpp(A, C, R,
  #                          kf$xitt, kf$xittm, kf$Ptt, kf$Pttm)

    xsmooth <- ks$xitT
    Vsmooth <- ks$PtT
    Wsmooth <- ks$PtTm

    delta <- matrix(0, os, ss)
    gamma <- matrix(0, ss, ss)
    beta <- matrix(0, ss, ss)

    for (t in 1:T) {
      z <- y[,t]; z[is.na(z)] <- 0
      # There could be an issue here
      delta <- delta + z %*% t(xsmooth[,t])
      gamma <- gamma + xsmooth[,t] %*% t(xsmooth[,t]) + Vsmooth[,,t]
      if (t > 1) {
        beta <- beta + xsmooth[,t] %*% t(xsmooth[,(t-1)]) + Wsmooth[,,t]
      }
    }

    gamma1 <- gamma - xsmooth[, T] %*% t(xsmooth[, T]) - Vsmooth[, , T]
    gamma2 <- gamma - xsmooth[, 1] %*% t(xsmooth[, 1]) - Vsmooth[, , 1]
    x1 <- xsmooth[, 1]
    V1 <- Vsmooth[, , 1]

    return(list(beta_t = beta, gamma_t = gamma, delta_t = delta, gamma1_t = gamma1,
                gamma2_t = gamma2, x1 = x1, V1 = V1, loglik_t = kf$loglik, xsmooth = xsmooth))

}



#' Implements a Kalman for dynamic factor model.
#'
#' @param initx Initial value for state space observations
#' @param initV Initial value for state covariance
#' @param x Observation matrix
#' @param A State space matrix
#' @param C System matrix
#' @param R State space covariance
#' @param Q System covariance
#' @return Filtered state space variable and its covariance matrix
#' as well as their forecast for next period for further iterations
#' @importFrom stats na.omit
K_filter <- function(initx, initV, x, A, C, R, Q) {
    T <- dim(x)[1]
    N <- dim(x)[2]
    r <- dim(A)[1]
    W <- !is.na(x)
    y <- t(x)

    xittm <- matrix(0, r, (T+1))
    xitt <- matrix(0, r, T)

    Pttm <- array(0, c(r, r, (T+1)))
    Ptt <- array(0, c(r, r, T))

    xittm[,1] <- initx
    Pttm[,,1] <- initV

    logl <- c()
    Ci <- C
    Ri <- R
    for (j in 1:T) {
#      missing_data <- MissData(y[,j], C, R)
#      C <- missing_data$C
#      R <- missing_data$R
      C <- Ci[W[j,],, drop=FALSE]
      R <- Ri[W[j,], W[j,], drop=FALSE]
      if (FALSE) #(all(!W[j,])) #(all(is.na(missing_data$y) == TRUE))
      {
         xitt[,,j] <- A %*% xittm[,,j]
         Ptt[,,j] <- C %*% Pttm[,,j] %*% t(C) + R
      } else
      {
         # Innovation covariance (inverse)
         Icov <- C %*% Pttm[,,j] %*% t(C) + R
         L <- solve(Icov)
         # Innovation residual
         Ires <- as.numeric(na.omit(y[,j])) - C %*% xittm[,j]
         # Optimal Kalman gain
         G <- Pttm[,,j] %*% t(C) %*% L
         # Updated state estimate: predicted + (Kalman gain)*fitted
         xitt[,j] <- xittm[,j] + G %*% Ires
         # Updated covariance estimate
         Ptt[,,j] <- Pttm[,,j] - G %*% C %*% Pttm[,,j]
         # State space variable and covariance predictions E[f_t | t-1]
         xittm[,(j+1)] <- A %*% xitt[,j]
         Pttm[,,(j+1)] <- A %*% Ptt[,,j] %*% t(A) + Q

         # Compute log-likelihood with Mahalanobis distance
         d <- length(Ires)
         S <- C %*% Pttm[,,j] %*% t(C) + R
         Sinv <- solve(S)
         if (nrow(R) == 1)
         {
           GG <- t(C) %*% solve(R) %*% C
           detS <- prod(R) %*% det(diag(1, r) + Pttm[,,j] %*% GG)
         } else {
           GG <- t(C) %*% diag(1/diag(R)) %*% C
           detS <- prod(diag(R)) * det(diag(1, r) + Pttm[,,j] %*% GG)
         }
         denom <- (2 * pi)^(d/2) * sqrt(abs(detS))
         mahal <- sum(t(Ires) %*% Ires) # Sinv %*% Ires)
         logl[j] <- -0.5 * mahal - log(denom)
        }
    }
    loglik <- sum(logl, na.rm=TRUE)
    return(list(xitt = xitt, xittm = xittm, Ptt = Ptt, Pttm = Pttm, loglik = loglik))
}

#' Implements Kalman smoothing and is used along with Kalman filter.
#' Kalman filter outputs enter Kalman smoother as inputs.
#'
#' @param A State space matrix
#' @param xitt State space variable
#' @param xittm Predicted state space variable
#' @param Ptt State space covariance
#' @param Pttm Predicted state space covariance
#' @param C System matrix
#' @param R State space covariance
#' @param W Logical matrix (T x n) indicating missing data.
#' TRUE if observation is present, FALSE if it is missing.
#' @return Smoothed state space variable and state space covariance matrix
K_smoother <- function(A, xitt, xittm, Ptt, Pttm, C, R, W) {
    T <- dim(xitt)[2]
    r <- dim(A)[1]

    Pttm <- Pttm[,,(1:(dim(Pttm)[3] - 1)), drop = FALSE]
    xittm <- xittm[,(1:(dim(xittm)[2] - 1)), drop = FALSE]

    # Whereas J is of constant dimension, L and K dimensions may vary
    # depending on existence of NAs
    J <- array(0, c(r, r, T))
    L <- list()
    K <- list()

    for (i in 1:(T-1)) {
        J[,,i] <- Ptt[,,i] %*% t(A) %*% solve(Pttm[,,(i+1)], tol = 1e-32)
    }

    Ci <- C
    Ri <- R
    for (i in 1:T) {
      # Only keep entries for non-missing data
      C <- Ci[W[i,],, drop=FALSE]
      R <- Ri[W[i,], W[i,], drop=FALSE]
      L[[i]] <- solve(C %*% Pttm[,,i] %*% t(C) + R)
      K[[i]] <- Pttm[,,i] %*% t(C) %*% L[[i]]
    }

    xitT <- cbind(matrix(0, r, (T-1)), xitt[,T])
    PtT <- array(0, c(r, r, T))
    PtTm <- array(0, c(r, r, T))
    PtT[,,T] <- Ptt[,,T]
    PtTm[,,T] <- (diag(1, r) - K[[T]] %*% C) %*% A %*% Ptt[,,(T-1)]

    for (j in 1:(T-1)) {
        xitT[,(T-j)] <- xitt[,(T-j)] + J[,,(T-j)] %*% (xitT[,(T+1-j)] - xittm[,(T+1-j)])
        PtT[,,(T-j)] <- Ptt[,,(T-j)] + J[,,(T-j)] %*% (PtT[,,(T+1-j)] - Pttm[,,(T+1-j)]) %*% t(J[,,(T-j)])
    }

    for (j in 1:(T-2)) {
        PtTm[,,(T-j)] <- Ptt[,,(T-j)] %*% t(J[,,(T-j-1)]) + J[,,(T-j)] %*% (PtTm[,,(T-j+1)] - A %*% Ptt[,,(T-j)]) %*% t(J[,,(T-j-1)])
    }

    return(list(xitT = xitT, PtT = PtT, PtTm = PtTm))
}


#' Implementation of Kim (1994) filter, an extension to Kalman filter
#' for dynamic linear models with Markov-switching. Documentation
#' is incomplete, rudimentary and needs to be rechecked!
#'
#' @param x0 Initial condition for state vector
#' @param P0 Initial condition for state variance
#' @param y Data matrix (Txn)
#' @param F System matrix for measurement equation
#' @param A Transition matrix for state equation
#' @param R Error covariance for measurement equation
#' @param Q Error covariance for state equation
#' @param p Transition probability matrix
#' @return Filtered states and covariances with associated probability matrices.
KimFilter <- function(x0, P0, y, F, A, R, Q, p)
{
  ## Define all containers for further computations. Notations for variables and
  ## indices, where appropriate, carefully follow Kim (1994). State vector is
  ## denoted as 'x', its covariance as 'P'. Appended letters explicit whether
  ## these are updated, approximated or smoothed.

  if (is.vector(y) || (ncol(y) <= length(x0)))
    stop("Number of factors should be strictly lower than number of variables. \n
          Increase number of variables or estimate a VAR model instead.")

  T <- nrow(y)
  n <- dim(F)[1]
  J <- length(x0)
  s <- dim(p)[1]

  ## x:   x^(i,j)_(t|t-1): predicted state vector - (2.6)
  ## xU:  x^(i,j)_(t|t): updated state vector - (2.11)
  ## P:   P^(i,j)_(t|t-1): predicted state covariance - (2.7)
  ## Pu:  P^(i,j)_(t|t): updated state covariance - (2.12)
  ## eta: eta^(i,j)_(t|t-1): conditional forecast error - (2.8)
  ## H:   H^(i,j)_(t): conditional variance of forecast error - (2.9)
  ## K:   K^(i,j)_(t): Kalman gain - (2.10)
  ## lik: f(y_t, S_(t-1)=i, S_t = j | t-1): joint conditional density - (2.16)
  ## loglik: log of (2.16)
  x <- array(NA, c(T,J,s,s))
  xU <- array(NA, c(T,J,s,s))
  P <- array(NA, c(T,J,J,s,s))
  Pu <- array(NA, c(T,J,J,s,s))
  eta <- array(NA, c(T,n,s,s))
  H <- array(NA, c(T,n,n,s,s))
  K <- array(NA, c(T,J,n,s,s))
  lik <- array(NA, c(T,s,s))
  loglik <- array(NA, c(T,s,s))
  ## Pr[S_(t-1) = i, S_t = j | t-1 ]: (2.15)
  ## Pr[S_(t-1) = i, S_t = j | t ]: (2.17)
  ## Pr[S_t = j | t-1 ]: used only for the smoothing part
  ## Pr[S_t = j | t ]: (2.18)
  jointP_fut <- array(NA, c(T,s,s))
  jointP_cur <- array(NA, c((T+1),s,s))
  stateP_fut <- array(NA, c(T,s))
  stateP <- array(NA, c(T,s))

  ## x^(j)_(t|t): approximate state vector conditional on S_j - (2.13)
  ## P^(j)_(t|t): approximate state covariance conditional on S_j - (2.14)
  xA <- array(NA, c(T,J,s))
  Pa <- array(0, c(T,J,J,s))
  result <- array(0, c(T,1))

  loglikf <- c()
  ## Some initial conditions to get started
  for (i in 1:s) { xA[1,,i] <- x0 }
  for (i in 1:s) { Pa[1,,,i] <- P0 }
  jointP_cur[1,,] <- matrix(c(0.25,0.25,0.25,0.25), ncol=2)
  temp <- array(NA, c(T,s,s))
  for (t in 2:T)
  {
    for (j in 1:s)
    {
      for (i in 1:s)
      {
        x[t,,i,j] <- A[,,j] %*% xA[(t-1),,i]
        P[t,,,i,j] <- A[,,j] %*% Pa[(t-1),,,i] %*% t(A[,,j]) + Q
        eta[t,,i,j] <- y[t,] - as(F[,,j], "matrix") %*% x[t,,i,j]
        H[t,,,i,j] <- F[,,j] %*% as(P[t,,,i,j], "matrix") %*% t(F[,,j]) + R
        K[t,,,i,j] <- P[t,,,i,j] %*% t(F[,,j]) %*% solve(H[t,,,i,j])
        xU[t,,i,j] <- x[t,,i,j] + K[t,,,i,j] %*% eta[t,,i,j]
        Pu[t,,,i,j] <- (diag(1,J) - K[t,,,i,j] %*% F[,,j]) %*% P[t,,,i,j]
        jointP_fut[t,i,j] <- p[i,j]*sum(jointP_cur[(t-1),,i]) # is everything alright here?
        lik[t,i,j] <- (2*pi)^(-n/2) * det(H[t,,,i,j])^(-1/2) *
                      exp(-1/2*t(eta[t,,i,j]) %*% solve(H[t,,,i,j]) %*% eta[t,,i,j]) *
                      jointP_fut[t,i,j]
        loglik[t,i,j] <- log(lik[t,i,j])
        jointP_cur[t,i,j] <- lik[t,i,j]
        loglikf[t] <- sum(loglik[t,,])
      }
      ## Technically, there should be sum(lik[t,,]) term but it cancels out and is computed later
      stateP[t,j] <- sum(jointP_cur[t,,j])
      stateP_fut[t,j] <- sum(jointP_fut[t,,j])
      ## Compute probability-filtered state process and its covariance
      xA[t,,j] <- xU[t,,,j] %*% jointP_cur[t,,j] / stateP[t,j]
      for (i in 1:s)
      {
        Pa[t,,,j] <- Pa[t,,,j] +
                     (Pu[t,,,i,j] + (xA[t,,j] - xU[t,,i,j]) %*% t(xA[t,,j] - xU[t,,i,j])) *
                     exp(log(jointP_cur[t,i,j]) - log(stateP[t,j]))
      }
    }
    jointP_cur[t,,] <- exp(log(jointP_cur[t,,]) - log(sum(lik[t,,])))
    stateP[t,] <- exp(log(stateP[t,]) - log(sum(lik[t,,])))
    result[t,1] <- log(sum(lik[t,,]))
  }

  return(list("result"=sum(result), "xA"=xA, "Pa"=Pa, "x"=x, "P"=P, "stateP"=stateP, "stateP_fut"=stateP_fut, "loglik"=loglik, "jointP"=jointP_fut, "lik"=result, "temp"=temp))

}

#' Smoothing algorithm from Kim (1994) to be used following a run
#' of KimFilter function.
#'
#' @param xA Filtered state vector to be smoothed
#' @param Pa Filtered state covariance to be smoothed
#' @param x State-dependent state vector
#' @param P State-dependent state covariance
#' @param A Array with transition matrices
#' @param p Markov transition matrix
#' @param stateP Evolving current probability matrix
#' @param stateP_fut Predicted probability matrix
#' @return Smoothed states and covariance matrices. This is the equivalent
#' of Kalman smoother in Markov-switching case.
KimSmoother2 <- function(xA, Pa, A, P, x, p, stateP, stateP_fut)
{
  ## Define all containers for further computations. Notations for variables and
  ## indices, where appropriate, carefully follow Kim (1994). State vector is
  ## denoted as 'x', its covariance as 'P'. Appended letters explicit whether
  ## these are updated, approximated or smoothed.

  T <- dim(xA)[1]
  J <- dim(xA)[2]
  s <- dim(xA)[3]

  ## Pr[S_t = j, S_(t+1) = k | T]: (2.20)
  ## Pr[S_t = j | T]: (2.21)
  jointPs <- array(NA, c(T,s,s))
  ProbS <- array(NA, c(T,s))

  ## xS: x^(j,k)_(t|T): inference of x_t based on full sample - (2.24)
  ## Ps: P^(j,k)_(t|T): covariance matrix of x^(j,k)_(t|T) - (2.25)
  ## Ptilde: helper matrix as defined after (2.25)
  xS <- array(0, c(T,J,s,s))
  Ps <- array(0, c(T,J,J,s,s))
  Ptilde <- array(NA, c(T,J,J,s,s))

  ## xAS: x^(j)_(t|T): smoothed and approximated state vector conditional on S_j (2.26)
  ## Pas: P^(j)_(t|T): smoothed and approximated state covariance conditional on S_j (2.27)
  ## xF: x_(t|T): state-weighted [F]inal state vector (2.28)
  ## Pf: P_(t|T): state-weighted [f]inal state covariance
  xAS <- array(0, c(T,J,s))
  Pas <- array(0, c(T,J,J,s))
  xF <- array(0, c(T,J))
  Pf <- array(0, c(T,J,J))
  # Initial conditions for smoothing loop
  ProbS[T,] <- stateP[T,]

  for (t in seq(T-1,1,-1))
  {
    for (j in 1:s)
    {
      for (k in 1:s)
      {
        jointPs[t,j,k] <- ProbS[(t+1),k]*stateP[t,j]*p[j,k] / stateP_fut[(t+1),k]
        Ptilde[t,,,j,k] <- Pa[t,,,j] %*% t(A[,,k]) %*% solve(P[(t+1),,,j,k])
        xS[t,,j,k] <- xA[t,,j] + Ptilde[t,,,j,k] %*% (xA[(t+1),,k] - x[(t+1),,j,k])
        Ps[t,,,j,k] <- Pa[t,,,j] +
                       Ptilde[t,,,j,k] %*% (Pa[(t+1),,,k] - P[(t+1),,,j,k]) %*% t(Ptilde[t,,,j,k])
        xAS[t,,j] <- xAS[t,,j] + jointPs[t,j,k]*xS[t,,j,k]
        Pas[t,,,j] <- Pas[t,,,j] + jointPs[t,j,k]*(Ps[t,,,j,k] +
                                          (xAS[t,,j] - xS[t,,j,k]) %*% t(xAS[t,,j] - xS[t,,j,k]))
      }
      ProbS[t,j] <- sum(jointPs[t,j,])
      xAS[t,,j] <- xAS[t,,j] / ProbS[t,j]
      Pas[t,,,j] <- Pas[t,,,j] / ProbS[t,j]
    }
  }
  for (t in 1:T)
  {
    for (j in 1:s)
    {
      xF[t,] <- xF[t,] + xAS[t,,j]*ProbS[t,j]
      Pf[t,,] <- Pf[t,,] + Pas[t,,,j]*ProbS[t,j]
    }
  }

  return(list("xF"=xF, "Pf"=Pf, "ProbS"=ProbS))
}
