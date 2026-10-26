#' Daily event and event-free counts
#'
#' Converts individual survival data into the daily sufficient statistics of
#' the Bayesian Kaplan-Meier model. Times are rounded up to integer days. An
#' event on day \eqn{t} takes priority over censoring on the same day, and a
#' participant censored on day \eqn{t} is counted as event-free through day
#' \eqn{t}.
#'
#' @param time Positive follow-up times (in days; non-integers are rounded up).
#' @param status Event indicator (1 = event, 0 = censored).
#' @param tau RMST time point \eqn{\tau} in days (positive integer). Only days
#'   \eqn{1,\dots,\tau} enter the RMST to \eqn{\tau}, so counts are returned
#'   for these days. \eqn{\tau} is distinct from the trial follow-up, which is
#'   not restricted. Use the same \eqn{\tau} for the historical and current
#'   trials.
#' @return An object of class \code{"bkm_counts"}: a list with \code{X}
#'   (number of events on each day \eqn{t = 1,\dots,\tau}), \code{Y} (number
#'   event-free, \eqn{N - X}), \code{N} (number at risk at the start of the
#'   day) and \code{tau}.
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
  X <- tabulate(t[ev], nbins = tau)
  structure(list(X = X, Y = N - X, N = N, tau = tau), class = "bkm_counts")
}

#' @export
print.bkm_counts <- function(x, ...) {
  cat("Daily counts over days 1 to", x$tau, ":", sum(x$X), "events,",
      x$N[1], "at risk on day 1\n")
  invisible(x)
}

as_counts <- function(x, tau) {
  if (is.null(x)) return(NULL)
  if (inherits(x, "bkm_counts")) {
    if (x$tau != tau) stop("counts object has a different RMST time point from 'tau'")
    return(x)
  }
  if (is.list(x) && all(c("time", "status") %in% names(x)))
    return(bkm_counts(x$time, x$status, tau))
  stop("each arm must be a 'bkm_counts' object or a data frame with 'time' and 'status'")
}

## internal: apply a function over replicates, in parallel on any platform
par_apply <- function(X, FUN, cores = 1) {
  if (cores <= 1) return(lapply(X, FUN))
  if (.Platform$OS.type == "windows") {
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(bkmborrow)))
    return(parallel::parLapply(cl, X, FUN))
  }
  parallel::mclapply(X, FUN, mc.cores = cores)
}
