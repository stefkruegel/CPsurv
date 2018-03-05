#' Simulate Survival Data with Change Point
#'
#' Simulates Weibull distributed survival data from a given data set with
#' change point above which hazard rate is constant.
#'
#' @inheritParams cpest
#' @param changeP Change point.
#' @param shape Shape parameter of Weibull distribution.
#' @param scale Scale parameter of Weibull distribution.
#' @param censoring Logical; if \code{TRUE}, censored data are generated.
#' @param censpoint Censoring point for Type I censoring.
#' @param times.int Logical; if \code{TRUE}, returned survival times are
#'   integers.
#' @param parametric Logical; if \code{TRUE}, survival times are generated
#'   parametrically by inverse transform sampling; otherwise Kaplan-Meier is
#'   used for simulation.

#' @return A dataset with survival times and corresponding censoring status
#'   ('event').
sim.survdata <-function(time, event, changeP, shape, scale, censoring, censpoint,
                        times.int, parametric){
  nobs <- length(time)

  if(parametric){
    # rate of exponential distr.
    rateE <- shape / scale * (changeP / scale) ^ (shape - 1)

    # Survivorfunction at change point
    St <- pweibull(q = changeP, shape = shape, scale = scale, lower.tail = FALSE)

    u <- runif(nobs)
    times.sim <- numeric(nobs)
    times.sim[u > St] <- qweibull(p = u[u > St], shape = shape, scale = scale,
                                  lower.tail = FALSE)
    times.sim[u <= St] <- changeP + qexp(p = u[u <= St] /
                                           pweibull(q = changeP, shape = shape,
                                                    scale = scale,
                                                    lower.tail = FALSE),
                                         rate = rateE, lower.tail = FALSE)
  } else{
    # nonparametric simulation of survivaltimes (by Kaplan-Meier)
    times.sim <- km.sim.survtimes(nobs = nobs, time = time, event = event,
                                  weibexp = TRUE, changeP = changeP)
  }

  data <- data.frame(time = times.sim, event = rep(1, nobs))

  #-----------------censoring--------------------
  if(censoring == "random"){
    censtimes <- time[event == 0]
    times.sim.cens <- km.sim.survtimes(nobs = nobs, time = censtimes,
                                       event = rep(1, length(censtimes)),
                                       weibexp = FALSE)

    data$event <- as.numeric(data$time <= times.sim.cens)
    data$time <- pmin(data$time, times.sim.cens)
  }else if(censoring == "type1"){
    data$event[data$time >= censpoint] <- 0
    data$time[data$time >= censpoint] <- censpoint
  }

  if(times.int){
    data$time <- ceiling(data$time)
  }
  data
}
