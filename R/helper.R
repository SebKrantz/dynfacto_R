#' Convergence test for EM-algorithm.
#'
#' @param loglik Current value of the log-likelihood function
#' @param previous_loglik Value of the log-likelihood function at the previous
#  iteration
#' @param threshold If difference is less than threshold, then algorithm has
#' converged
#' @param check_increased TO DOCUMENT
#' @return A logical statement indicating whether EM algorithm has converged
#' according to slope convergence test
em_converged <- function(loglik, previous_loglik, threshold=1e-4, check_increased=TRUE) {

    converged <- FALSE
    decrease <- 0

    if (check_increased == TRUE) {
        if (loglik - previous_loglik < -0.001) {
#            message("*** Likelihood decreased from ", previous_loglik, " to ", loglik, "\n")
            decrease <- 1
        }
    }
    if (loglik != Inf & previous_loglik != Inf) {
    delta_loglik <- abs(loglik - previous_loglik)
    avg_loglik <- (abs(loglik) + abs(previous_loglik) + .Machine$double.eps)/2

    if ((delta_loglik/avg_loglik) < threshold) {
        converged <- TRUE
    }
    return(converged) }
    else return(FALSE)
    # return(list('converged'=converged, 'decrease'=decrease))
}


#' Estimate a p-th order vector autoregressive (VAR) model
#'
#' @param x Data matrix (T x n)
#' @param p Maximum lag order, i.e. VAR(p) will be estimated
#' @return Estimated parameter matrix, residuals and regression model
#' independent and dependent variables
#' @export
#' @examples
#' x <- matrix(rnorm(50*2), nrow=50, ncol=2)
#' VAR(x, 2)
VAR <- function(x, p) {
    T <- nrow(x)
    Y <- x[(p + 1):T, ]
    X <- c()
    for (i in 1:p) {
        X <- cbind(X, x[(p + 1 - i):(T - i), ])
    }
    A <- solve(t(X) %*% X) %*% t(X) %*% Y
    res <- Y - X %*% A

    return(list(Y = Y, X = X, A = A, res = res))
}


generateArray <- function(x) {

  dims1 <- dim(x)
  dims2 <- dim(x[1,1][[1]])

  nDims <- length(dims1) + length(dims2)
  array(unlist(x), c(dims2, dims1))

}
