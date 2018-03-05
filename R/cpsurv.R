#' @title Nonparametric Change Point Estimation
#'
#' @description Change point estimation for survival data based on exact
#'   binomial test.
#'
#' @details Change point is a point in time, from which on the hazard rate is
#'   supposed to be constant. For its estimation the timeline up to \code{cpmax}
#'   is split into equidistant intervals of width \code{intwd} and exact
#'   binomial tests are executed for each interval. The change point is
#'   estimated by fitting a regression model on the resulting p-values. See
#'   Brazzale \emph{et al} (2017) for details. \cr\cr For bootstrap bias
#'   correction the change point is estimated for a given number
#'   (\code{B.correct}) of bootstrap samples whereupon the bias is built by
#'   subtracting their median from primary estimation. Depending on argument
#'   \code{parametric} the data for bootstrapping are simulated either
#'   parametric (Weibull distributed with estimated shape and scale parameters)
#'   or nonparametric (based on Kaplan-Meier estimation).
#' @inheritParams cpest
#' @param censoring Type of right-censoring for simulated data on which the
#'   bootstrap bias correction is based. Possible types are "random" for
#'   \emph{random censoring} (default), "type1" for \emph{Type I censoring} or
#'   "no" for data without censored observations. Because simulated data should
#'   be similar to given data, the censoring type is adapted from vector
#'   'events' if given and argument 'censoring' is ignored than.
#' @param censpoint Point of \emph{Type I censoring}; if missing, minimum time
#'   after which all events are equal to 0 is used. Censpoint is only needed for
#'   bootstrap bias correction.
#' @param biascorrect Logical; if \code{TRUE}, a bootstrap bias correction is
#'   performed; see 'Details'.
#' @param parametric Indicator for parametric bias-correction (see Details for
#'   more information).
#' @param B.correct Number of bootstrap samples for bias-correction; defaults to
#'   49.
#' @param opt.start Numeric vector of length two; initial values for the Weibull
#'   parameters (shape and scale parameters) to be optimized if parametric
#'   bootstrap bias correction is used.
#' @param boot.ci Indicator if confidence intervals (and thereby standard
#'   deviation) should be calculated by bootstrap sampling. Please note the
#'   extended runtime (see details for examples).
#' @param B Number of bootstrap samples for confidence intervals; defaults to
#'   999.
#' @param conf.level Confidence level for bootstrap confidence intervals.
#' @param seed Seed for random number generator (optional).
#' @param parallel Indicator if bootstrap-sampling is executed parallelized
#'   (based on package 'parallel'); operating system is identified
#'   automatically.
#' @param cores Number of CPU-cores that are used for parallelization; maximum
#'   possible value is the detected number of logical CPU cores.
#' @return \tabular{ll}{ \code{cp}\tab estimated change point\cr
#'   \code{p.values}\tab p-values resulting from exact binomial test\cr
#'   \code{pv.mean}\tab mean of p-values for intervals above the estimated
#'   change point\cr \code{lower.lim}\tab lower interval limits\cr
#'   \code{upper.lim}\tab upper interval limits\cr \code{cp.bc}\tab bias
#'   corrected change point\cr \code{ml.shape}\tab ML estimator of shape
#'   parameter for Weibull distribution\cr \code{ml.scale}\tab ML estimator of
#'   scale parameter for Weibull distribution\cr \code{cp.boot}\tab estimated
#'   change points for bootstrap samples\cr \code{sd}\tab standard deviation
#'   estimated by bootstrap sampling\cr \code{ci.normal}\tab confidence interval
#'   with normal approximation\cr \code{ci.percent}\tab bootstrap percentile
#'   interval\cr \code{conf.level}\tab the \code{conf.level} argument passed to
#'   \code{cpsurv}\cr \code{B}\tab the \code{B} argument passed to
#'   \code{cpsurv}\cr \code{time}\tab the \code{time} argument passed to
#'   \code{cpsurv}\cr \code{event}\tab the \code{event} argument passed to
#'   \code{cpsurv}\cr \code{cpmax}\tab the \code{cpmax} argument passed to
#'   \code{cpsurv}\cr \code{intwd}\tab the \code{intwd} argument passed to
#'   \code{cpsurv}\cr \code{call}\tab matched call}
#' @references Brazzale, A. R. and Küchenhoff, H. and Krügel, S. and Hartl, W.
#'   (2017) \emph{Nonparametric change point estimation for survival
#'   distributions with a partially constant hazard rate.}
#' @author Stefanie Krügel \email{stefanie.kruegel@@gmail.com}
#' @examples
#' data(survdata)
#' # estimate change point for survdata (random censored)
#' cp <- cpsurv(survdata$time, survdata$event, cpmax = 360, intwd = 20)
#' summary(cp)
#'
#' \dontrun{
#' # estimation with parametric bootstrap bias correction
#' cp_param <- cpsurv(survdata$time, survdata$event, cpmax = 360, intwd = 20,
#'             biascorrect = TRUE, parametric = TRUE)
#' summary(cp_param)
#'
#' # with bootstrap confidence intervals and parametric bootstrap bias
#' cp_ci <- cpsurv(survdata$time, survdata$event, cpmax = 360, intwd = 20,
#' biascorrect = TRUE, parametric = FALSE, boot.ci = TRUE, cores = 4, seed = 36020)
#' # runtime: approx. 180 min (with Intel(R) Core(TM) i7 CPU 950 @ 3.07GHz, 4 logical CPUs used)
#' }
#' @import parallel
#' @importFrom stats dweibull median nlminb pbinom pexp predict pweibull qexp
#'   qnorm qweibull runif sd smooth.spline spline
#' @export
cpsurv <- function(time, event, cpmax, intwd, cpmin = 0,
                   censoring = c("random", "type1", "no"), censpoint = NULL,
                   biascorrect = FALSE, parametric = FALSE, B.correct = 49,
                   opt.start = c(0.1, 50), boot.ci = FALSE,
                   B = 999, conf.level = 0.95, norm.riskset = TRUE, seed = NULL,
                   parallel = TRUE, cores = 4L){

  if (!is.null(seed)) {
    if(parallel) RNGkind(kind = "L'Ecuyer-CMRG")
    set.seed(as.integer(seed))
  }

  # function for input checking
  is_single_num <- function(x, lower = -Inf, upper = Inf){
    if(is.na(x) || !is.numeric(x) || (length(x) != 1)){
      stop("Argument '", substitute(x), "' is not a single numeric value.")
    }
    if(x < lower) stop("Value for argument '", substitute(x), "' too low.")
    if(x > upper) stop("Value for argument '", substitute(x), "' too high.")
  }

  stopifnot(is.numeric(time),
            all(sapply(time, function(z) z >= 0)),
            length(time) >= 1)

  if(missing(event))
    event <- rep.int(1, length(time))
  stopifnot(is.numeric(event), length(event) >= 1)
  if(!all(sapply(event, function(z) z %in% c(0,1)))){
    stop("Argument 'event' has to be binary.")
  }

  is_single_num(cpmax, lower = 1L)
  is_single_num(cpmin, lower = 0L, upper = cpmax)
  if(!is.null(censpoint)) is_single_num(censpoint, lower = 0L, upper = max(time))
  is_single_num(B.correct, lower = 1L)
  is_single_num(B, lower = 1L)
  is_single_num(conf.level, lower = 0L, upper = 1L)
  is_single_num(cores, lower = 1L)
  if(missing(intwd)) intwd <- ceiling(cpmax / 20)
  is_single_num(intwd, lower = 1L)

  stopifnot(!is.na(norm.riskset), is.logical(norm.riskset), length(norm.riskset) == 1,
            !is.na(biascorrect), is.logical(biascorrect), length(biascorrect) == 1,
            !is.na(parametric), is.logical(parametric), length(parametric) == 1,
            !is.na(boot.ci), is.logical(boot.ci), length(boot.ci) == 1,
            !is.na(parallel), is.logical(parallel), length(parallel) == 1)

  if((cpmax%%intwd) != 0){
    cpmax <- cpmax + (cpmax %% intwd)
    warning("'cpmax' is not a multiple of 'intwd'; cpmax corrected to ", cpmax)
  }

  censoring <- match.arg(censoring)

  if(!identical(length(time), length(event)))
    stop("Vectors 'time' and 'event' must be of equal length.")

  if(sum(event[time > cpmax] == 1) == 0){
    stop("No events with 'time' > 'cpmax'; ",
         "choose lower value for 'cpmax'.")
  }

  if(!is.vector(opt.start) || !length(opt.start) == 2 || !all(opt.start > 0)){
    stop("Invalid 'opt.start' argument. See documentation.")
  }
  #----------------------------------------------------------------------------

  # check censoring
  if(all(event %in% 1)){
    censoring <- "no"
  }
  if(censoring == "type1" && is.null(censpoint)){
    # point of type 1 censoring
    censpoint <- min(time[sapply(seq_along(time), function(i)
      all(event[time >= time[i]] == 0))])
    warning("point for Type I censoring missing; 'censpoint' was set to ", censpoint)
  }

  nobs <- length(time)

  # if times of given data are integers, simulated times are also integers
  check.integer <- function(x){
    x == round(x)
  }
  if(all(sapply(time, check.integer))){
    times.int <- TRUE
  } else{
    times.int <- FALSE
  }

  # estimated change point
  cp <- cpest(time = time, event = event, intwd = intwd, cpmax = cpmax,
              cpmin = cpmin, norm.riskset = norm.riskset)
  changeP <- cp$cp

  if(biascorrect){
    bcresult <- bootbiascorrect(changeP = changeP, time = time, event = event,
                                censoring = censoring, censpoint = censpoint,
                                intwd = intwd, opt.start = opt.start,
                                cpmax = cpmax, cpmin = cpmin,
                                norm.riskset = norm.riskset,
                                B.correct = B.correct, parametric = parametric,
                                times.int = times.int)
    if(is.na(bcresult$bc)){
      if(parametric == TRUE){
        warning("Parametric bootstrap bias correction not possible;",
                "maybe times are not Weibull distributed")
      } else{
        warning("Nonparametric bootstrap bias correction not possible",
                "(no events above 'cpmax' in simulated data)")
      }
    }
    cp$cp.bc <- bcresult$bc
    if(parametric){
      cp$ml.shape <- bcresult$shape
      cp$ml.scale <- bcresult$scale
    }
    if(bcresult$oob[1] == TRUE){
      warning("bias corrected change point out of bounds; it was set to 'cpmin'")
    } else if(bcresult$oob[2] == TRUE){
      warning("bias corrected change point out of bounds; it was set to 'cpmax'")
    }
  }


  # variance estimation by bootstrap sampling
  if(boot.ci){

    bootstrap <- function(x){
      boot.out <- NA

      while(is.na(boot.out)){
        samp <- sample(1:nobs, size = nobs, replace=T)
        cp.tmp <- cpest(time = time[samp], event = event[samp], intwd = intwd,
                        cpmax = cpmax, cpmin = cpmin, norm.riskset = norm.riskset)$cp
        if(biascorrect){
          boot.out <- bootbiascorrect(changeP = cp.tmp, time = time[samp],
                                      event = event[samp], censoring = censoring,
                                      censpoint = censpoint, intwd = intwd,
                                      cpmax = cpmax, cpmin = cpmin,
                                      norm.riskset = norm.riskset, B.correct = B.correct,
                                      parametric = parametric, times.int = times.int,
                                      opt.start = opt.start)$bc
        } else{
          boot.out <- cp.tmp
        }
      }
      boot.out
    }

    if(parallel){
      if(parallel::detectCores() < cores) cores <- parallel::detectCores()
      if(.Platform$OS.type != "windows") {
        cpboot <- parallel::mclapply(1:B, FUN = bootstrap, mc.cores = cores)
      } else{
        cl <- parallel::makePSOCKcluster(rep("localhost", cores))
        parallel::clusterExport(cl = cl, varlist = c("bootbiascorrect",
          "neg.loglik.WeibExp", "cpest", "sim.survdata",
          "km.sim.survtimes", "nobs", "time", "event", "intwd", "cpmax",
          "cpmin", "norm.riskset", "censoring", "censpoint", "biascorrect", "B.correct",
        "parametric", "times.int", "opt.start"),
          envir = environment())

        if(!is.null(seed)) parallel::clusterSetRNGStream(cl)
        cpboot <- parallel::parLapply(cl, 1:B, fun = bootstrap)
        parallel::stopCluster(cl)
      }
      cpboot <- unlist(cpboot)
    } else{
      cpboot <- vapply(1:B, bootstrap, FUN.VALUE = numeric(1L))
    }

    cp$cp.boot <- cpboot

    # standard deviation
    sd_boot <- sd(cpboot)
    cp$sd <- sd_boot

    # confidence intervals
    quant <- (1 + c(-conf.level, conf.level))/2

    #normal approximation
    if(biascorrect){
      ci.normal <- c(cp$cp.bc - qnorm(quant[2]) * sd_boot,
                     cp$cp.bc + qnorm(quant[2]) * sd_boot)
    } else{
      ci.normal <- c(changeP - qnorm(quant[2]) * sd_boot,
                     changeP + qnorm(quant[2]) * sd_boot)
    }

    if(ci.normal[1] < 0){
      ci.normal[1] <- 0
      warning("lower bound of confidence interval was set to 0")
    }
    if(ci.normal[2] > cpmax){
      ci.normal[2] <- cpmax
      warning("upper bound of confidence interval was set to 'cpmax'")
    }
    cp$ci.normal <- ci.normal

    #percentile interval
    perc <- (B + 1) * quant
    ptrunc <- trunc(perc)
    cpsort <- sort(cpboot)
    percint <- numeric(2L)
    percint[1] <- ifelse(ptrunc[1] == 0, cpsort[1L], cpsort[ptrunc[1]])
    percint[2] <- ifelse(ptrunc[2] == B, cpsort[B], cpsort[ptrunc[2]])
    if(ptrunc[1] == 0 || ptrunc[2] == B){
      warning("extreme order statistics used as endpoints for percentile interval")
    }
    cp$ci.percent <- percint

    cp$conf.level <- conf.level
    cp$B <- B
  }

  structure(c(cp, list(time = time, event = event, cpmax = cpmax, intwd = intwd,
                       call = match.call())), class = "cpsurv")
}
