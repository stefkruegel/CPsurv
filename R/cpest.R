#' Estimates change point using shifted intervals
#'
#' Shifts intervals iteratively and estimates change point at each step. Final
#' change point is calculated by optimization over all estimations.
#'
#' @param time Numeric vector with survival times.
#' @param event Numeric vector indicating censoring status; 0 = alive (censored), 1 =
#'   dead (uncensored). If missing, all observations are assumed to be
#'   uncensored.
#' @param intwd Width of intervals into which the time period is split; default
#'   is \code{ceiling(cpmax/20)}. Has to be an integer value.
#' @param cpmax Upper bound for estimated change point. Time period is split into
#'   intervals up to this point. Has to be an integer value.
#' @param cpmin Lower bound for estimated change point; default is
#'   \code{cpmin=0}. Has to be an integer value.
#' @param norm.riskset Logical; if \code{TRUE} normalized number of units at
#'   risk is used within an interval.
#' @return A list with estimated change point, p-values of exact binomial
#'   test, mean of p-values above estimated change point (part of regression
#'   function), lower and upper bounds of confidence intervals.
#' @seealso \code{\link{cpsurv}}

cpest <- function(time, event, cpmax, intwd, cpmin, norm.riskset){
  nobs <- length(time)

  # estimation of lambda (const. hazard rate for t > cpmax)
  rate <- sum(event[time > cpmax] == 1) / sum(time[time > cpmax] - cpmax)

  #probability for event within an interval
  pr = pexp(intwd, rate)

  # determine interval limits
  limseq <- seq(cpmin, (cpmax + 2 * intwd - 1))
  limseq <- limseq[1:(intwd * floor(length(limseq) / intwd))]
  lim <- matrix(limseq, nrow = intwd)

  results <- matrix(nrow = intwd, ncol = 3)
  pvals <- matrix(nrow = intwd, ncol = ncol(lim) - 1)

  for(shift in 1:intwd){
    # number of events within intervals
    x <- table(cut(time[event == 1], lim[shift,]))

    if(shift == 1){
      cumtab <- c(0, cumsum(table(cut(time, lim[shift,]))))
    } else{
      cumtab <- cumsum(table(cut(time, c(0, lim[shift,]))))
    }
    n <- nobs - cumtab[-length(cumtab)]

    if(norm.riskset){
      # times of right-censored data
      tcens <- time[event == 0]
      censcut <- cut(tcens, lim[shift, ])
      splitted <- split(tcens, censcut)
      # correct n by censored observations
      n <- vapply(1:(ncol(lim) - 1), function(k){
        n[k] - round(sum((lim[shift, k + 1] - splitted[[k]]) / intwd))},
        FUN.VALUE = numeric(1L))
    }

    pv <- vapply(1:(ncol(lim) - 1), function(i){
      # pv = P(k >= x) = P(k > x-1)
      pbinom(x[i] - 1, n[i], pr, lower.tail = FALSE)

    }, FUN.VALUE = numeric(1L))

    # means of p-values
    mean_pv <- rev(cumsum(rev(pv))) / rev(seq_along(pv))
    lower <- lim[shift, -ncol(lim)]
    S <- vapply(1:length(lower),
                function(i){sum((pv - mean_pv[i] * (lower[i] <= lower))^2)},
                FUN.VALUE = numeric(1L))

    opt <- which.min(S)

    cp <- lim[shift, opt]
    if(cp < cpmin) cp <- cpmin
    if(cp > cpmax) cp <- cpmax

    results[shift, ] <- c(min(S), cp, mean_pv[opt])
    pvals[shift, ] <- pv
  }
  optrow <- which.min(results[ ,1])

  list(cp = results[optrow, 2],
       p.values = pvals[optrow, ],
       pv.mean = results[optrow, 3],
       lower.lim = lim[optrow, -ncol(lim)],
       upper.lim = lim[optrow, -1],
       rate = rate)
}
