#' Implements Bootstrap Bias Correction
#'
#' @inheritParams cpest
#' @param changeP Estimated change point.
#' @param censoring Type of right-censoring for simulated data on which the
#'   bootstrap bias correction is based. Possible types are "random" for
#'   \emph{random censoring} (default), "type1" for \emph{Type I censoring} or
#'   "no" for data without censored observations. Because simulated data should
#'   be similar to given data, the censoring type is adapted from vector
#'   'events' if given and argument 'censoring' is ignored than.
#' @param censpoint Point of \emph{Type I censoring}; if missing, minimum time
#'   after which all events are equal to 0 is used. Censpoint is only needed for
#'   bootstrap bias correction.
#' @param B.correct Number of bootstrap samples for bias correction; defaults to
#'   49.
#' @param parametric Logical; if \code{TRUE} parametric bootstrap bias
#'   correction is used (simulation of boostrap samples is based on estimated
#'   Weibull parameters); otherwise Kaplan-Meier is used for a nonparametric
#'   bootstrap bias correction.
#' @param times.int Logical; if \code{TRUE} simulated survival times are
#'   integers.
#' @param opt.start Numeric vector of length two; initial values for the Weibull
#'   parameters (shape and scale parameters) to be optimized if parametric bootstrap
#'   bias correction is used.
#' @return A list with bias-corrected change point and optional estimated shape
#'   and scale parameters of the Weibull distribution.
#' @importFrom survival Surv survfit
bootbiascorrect <- function(changeP, time, event, censoring, censpoint, intwd,
                            cpmax, cpmin, norm.riskset, B.correct, parametric,
                            times.int, opt.start){

  # for parametric bootstrap bias correction estimate shape and scale parameters
  # of Weibull distribution
  if(parametric){
    par.tmp <- nlminb(start = opt.start, objective = neg.loglik.WeibExp,
                      changeP = changeP, time = time, event = event,
                      lower = c(1e-08, 1e-08))$par
    shape <- par.tmp[1]
    scale <- par.tmp[2]
  }


  # simulate data and estimate change point 'B.correct' times
  cp.boot <- vapply(1:B.correct, function(idx){

    # simulate data with estimated parameters
    data.boot <- sim.survdata(time = time, event = event, changeP = changeP,
                              censoring = censoring, censpoint = censpoint,
                              scale = scale, shape = shape,
                              times.int = times.int, parametric = parametric)

    if(sum(data.boot[data.boot$time > cpmax, "event"] == 1) == 0){
      # NA is returned, if there are no events above 'cpmax'
      NA
    } else{
      # return estimated change point for current bootstrap-sample
      cpest(time = data.boot[,"time"], event = data.boot[,"event"],
            intwd = intwd, cpmax = cpmax, cpmin = cpmin,
            norm.riskset = norm.riskset)$cp
    }
  }, numeric(1L))

  # indicator, if bc is out of bounds
  oob <- logical(2L)

  if(any(is.na(cp.boot))){
    if(parametric){
      return(list(bc = NA, shape = shape, scale = scale, oob = oob))
    } else return(list(bc = NA, oob = oob))
  }

  bias <- median(cp.boot) - changeP
  bc <- changeP - bias


  if(bc < cpmin){
    bc <- cpmin
    oob[1] <- TRUE
  }
  if(bc > cpmax){
    bc <- cpmax
    oob[2] <- TRUE
  }

  if(parametric){
    list(bc = bc, shape = shape, scale = scale, oob = oob)
  } else{
    list(bc = bc, oob = oob)
  }
}
