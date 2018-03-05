#' @name survdata
#' @title Simulated Survival Data
#' @description A simulated dataset with 1500 fake right-censored survival times
#'   with a change point at \code{time = 90}.\cr The survival times are Weibull
#'   distributed with parameters \code{shape = 0.44} and \code{scale = 100}
#'   below the change point and have a constant hazard rate above.
#' @docType data
#' @usage survdata
#' @format \tabular{ll}{ \code{time} \tab survival or censoring time\cr
#' \code{event} \tab censoring status (0 = alive, 1 = dead)}
NULL
