#' Negative Log-Likelihood for Weibull-Exponential Distribution
#'
#' @param param Shape and scale parameter for Weibull distribution.
#' @param changeP Changepoint.
#' @param time Vector of survival times.
#' @param event Vector indicating censoring status; 0 = alive (censored), 1 = dead
#'   (uncensored).
#' @return Value of the negative log-likelihood.
neg.loglik.WeibExp <- function(param, changeP, time, event)
{
  shape <- param[1]
  scale <- param[2]
  rateE <- shape/scale * (changeP/scale)^(shape-1)
  #
  is.W <- as.numeric(time <= changeP)
  x.E <- (time - changeP)
  #
  ll1 <- event * is.W * dweibull(x=time, shape=shape, scale=scale)
  ll2 <- (1-event) * is.W * pweibull(q=time, shape=shape, scale=scale, lower.tail=F)
  ll3 <- event * (1-is.W) * rateE * pweibull(q=changeP, shape=shape, scale=scale,
                                           lower.tail=F) * pexp(q=x.E, rate=rateE, lower.tail=F)
  ll4 <- (1-event) * (1-is.W) * pweibull(q=changeP, shape=shape, scale=scale,
                                       lower.tail=F) * pexp(q=x.E, rate=rateE, lower.tail=F)
  ll <- ll1 + ll2 + ll3 + ll4
  -sum(log(ll))
}
