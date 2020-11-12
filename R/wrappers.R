#' Kim filtering
#'
#' @return This function returns the same output as \code{\link{KimFilterCpp}}.
#'   However, it converts certain elements of the Kim Filter to standard multi-dimensional
#'   R arrays and so its outputs can be used as inputs for \code{\link{KimSmoother}}.
#' @inheritParams KimFilterCpp
#' @seealso \code{\link{KimFilterCpp}}
#' @export
KimFilter <- function(y, R, Q, F1, A1, x0, P0, p) {
  kf <- KimFilterCpp(y, R, Q, F1, A1, x0, P0, p)

  kf$x <- generateArray(kf$x)
  kf$P <- generateArray(kf$P)
  kf$Pa <- generateArray(kf$Pa)

  kf
}
