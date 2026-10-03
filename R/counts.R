#' Daily event and event-free counts
#'
#' Converts individual survival data into the daily sufficient statistics of
#' the discrete-time (daily) Bayesian Kaplan-Meier model. Times are rounded up
#' to integer days. An event on day \eqn{d} takes priority over censoring on
#' the same day, and a participant censored on day \eqn{d} is counted as
#' event-free through day \eqn{d}.
#'
#' @param time Positive follow-up times (in days; non-integers are rounded up).
#' @param status Event indicator (1 = event, 0 = censored).
#' @param tau Restriction horizon in days (positive integer).
#' @return An object of class \code{"bkm_counts"}: a list with \code{D}
#'   (events on each day \eqn{1,\dots,\tau}), \code{Y} (at risk minus events),
#'   \code{N} (at risk) and \code{tau}.
#' @examples
#' bkm_counts(time = c(3, 5, 5, 8), status = c(1, 0, 1, 1), tau = 6)
#' @export
bkm_counts <- function(time, status, tau) {
  if (length(time) != length(status)) stop("'time' and 'status' differ in length")
  if (any(!is.finite(time)) || any(time <= 0)) stop("'time' must be positive and finite")
  if (!all(status %in% c(0, 1))) stop("'status' must be 0/1")
  tau <- as.integer(tau)
  if (length(tau) != 1L || tau < 1L) stop("'tau' must be a positive integer")
  t <- as.integer(ceiling(time))
  ev <- status == 1 & t <= tau
  x <- pmin(t, tau)
  N <- rev(cumsum(rev(tabulate(x, nbins = tau))))
  D <- tabulate(t[ev], nbins = tau)
  structure(list(D = D, Y = N - D, N = N, tau = tau), class = "bkm_counts")
}

#' @export
print.bkm_counts <- function(x, ...) {
  cat("Daily counts over", x$tau, "days:", sum(x$D), "events,",
      x$N[1], "at risk on day 1\n")
  invisible(x)
}

as_counts <- function(x, tau) {
  if (is.null(x)) return(NULL)
  if (inherits(x, "bkm_counts")) {
    if (x$tau != tau) stop("counts object has a different horizon from 'tau'")
    return(x)
  }
  if (is.list(x) && all(c("time", "status") %in% names(x)))
    return(bkm_counts(x$time, x$status, tau))
  stop("each arm must be a 'bkm_counts' object or a data frame with 'time' and 'status'")
}
