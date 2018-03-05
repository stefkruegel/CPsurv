#' @include cpsurv.R
NULL

#' @title Summarize and print cpsurv objects
#' @description Summary and print methods for objects inheriting from a call to
#'   \code{\link{cpsurv}}.
#' @seealso \code{\link{cpsurv}}
#' @name summarize.cpsurv
NULL

#' @param x An object of class \code{cpsurv} or \code{summary.cpsurv} to be
#'   printed out.
#' @param ... not used
#' @rdname summarize.cpsurv
#' @export
print.cpsurv <- function(x, ...){
  cat("Estimated change point:\n")
  if(!is.null(x$cp.bc)){
    cat(" ", x$cp.bc, "\n")
  }else{
    cat(" ", x$cp, "\n")
  }
}

#' @details The main results from \code{cpsurv} are printed out in a
#'   well-arranged format. If the estimated change point is bias corrected, both
#'   estimates (the original, and the corrected one) are shown in the summary.
#'   If a bootstrap-sampling was executed, the output contains a summary of the
#'   resultant bootstrap-estimates.
#' @param object An object of class \code{cpsurv}.
#' @examples
#' data(survdata)
#' cpest <- cpsurv(survdata$time, survdata$event, cpmax = 360)
#' summary(cpest)
#' @rdname summarize.cpsurv
#' @export
summary.cpsurv <- function(object, ...){
  if(is.null(object$cp.bc)){
    estimate <- object$cp
  } else{
    estimate <- object$cp.bc
  }

  if(!is.null(object$sd)){
    estimate["std.error"] <- round(object$sd, 3)

    qq <- stats::quantile(object$cp.boot)
    names(qq) <- c("Min", "1Q", "Median", "3Q", "Max")

    quant <- (1 + c(-object$conf.level, object$conf.level))/2

    normal <- as.character(round(object$ci.normal, 4))
    names(normal) <- c(paste0((quant[1]* 100), "%"), paste0((quant[2] * 100), "%"))
    percentile <- as.character(object$ci.percent)

    int.out <- rbind(normal, percentile)

    colnames(int.out) <- c(paste0((quant[1]* 100), "%"), paste0((quant[2] * 100), "%"))

    structure(c(object,
                list(estimate = estimate, qq = qq, int.out = int.out, normal = normal,
                     percentile = percentile)),
              class = "summary.cpsurv")
  } else{
    structure(c(object, list(estimate = estimate)),
              class = "summary.cpsurv")
  }
}

#' @rdname summarize.cpsurv
#' @export
print.summary.cpsurv <- function(x, ...){
  bc <- !is.null(x$cp.bc)
  boot <- !is.null(x$sd)
  param <- !is.null(x$ml.shape)

  cat("\nNONPARAMETRIC CHANGE POINT ESTIMATION\n")
  #cat("``````````````````´´´´´´´´´´´´´´´´´´\n")
  if(bc){
    partext <- ifelse(param, "PARAMETRIC", "NONPARAMETRIC")
    cat("- with", partext, "Bootstrap Bias Correction\n")
  }
  cat("\nCALL: \n")
  dput(x$call, control = NULL)


  if(boot){
    names(x$estimate) <- c("change point", "std.error")
    cat("\nEstimation: \n")
    print(x$estimate)
  } else{
    cat("\nestimated change point: ", x$estimate, "\n")
  }
  cat("(based on exact binomial tests for intervals of width ",
      x$intwd, ")\n", sep = "")

  if(boot){
    cat("\n\nBootstrap Confidence Intervals:\n")
    print(x$int.out, quote = FALSE)
    cat("\n(based on", x$B, "bootstrap replicates)\n")
    cat("\nBootstrap-Estimations:\n")
    print(x$qq)
  }
  if(param){
    cat("\n\nEstimated Weibull parameters (used for bias correction):\n")
    cat("   shape:", round(x$ml.shape, 3), "    scale:", round(x$ml.scale, 3))
  }
}


#' Plot method for objects of class cpsurv
#'
#' Plot method for objects of class 'cpsurv' inheriting from a call to
#' \code{\link{cpsurv}}.
#'
#' @details The value \code{type = "pvals"} produces a plot with p-values used
#'   to estimate the stump regression model with superimposed least squares
#'   regression line. For \code{type = "events"} a barplot is produced with
#'   frequency of events per unit at risk for each interval (with length
#'   \code{intwd}. For \code{type = "hazard"} the estimated hazard rate (based
#'   on \code{\link[muhaz]{muhaz}}) is plotted with optional (normal- or
#'   percentile-) confidence intervals and the estimated constant hazard rate.
#'
#' @param x An object of class 'cpsurv' (estimated with \code{cpsurv}).
#' @param type A vector of character strings to select the plots for printing.
#'   The value should be any subset of the values c("pvals", "events", "hazard")
#'   or simply "all", where all possible plots are shown.
#' @param ci.type Character representing the type of confidence interval to plot
#'   (if existing); "perc" for percentile interval and "norm" for CI with normal
#'   approximation (default is "perc").
#' @param ci Logical; if \code{TRUE}, a bootstrap confidence interval is plotted
#'   (if existing).
#' @param const.haz Logical; if \code{TRUE}, the estimated constant hazardrate
#'   is plotted.
#' @param regline Logical; if \code{TRUE}, the regression line is plotted.
#' @param legend Logical; if \code{TRUE}, the plots contain legends.
#' @param xlim Vector with x limits (timeline) for each plot if supplied;
#'   default is c(0, x$cpmax).
#' @param ylim Vector with y limits for plots of type "events" and "hazard". For
#'   changing ylim for only one of them, plot them separately by use of argument
#'   'type'.
#' @param main Main title for each plot if supplied.
#' @param xlab Character vector used as x label for all plots if supplied.
#' @param ylab Character vector used as y label for all plots if supplied.
#' @param min.time Left bound of time domain used for
#'   \code{\link[muhaz]{muhaz}}. If missing, min.time is considered 0.
#' @param max.time Right bound of time domain used for
#'   \code{\link[muhaz]{muhaz}}. If missing, value 'cpmax' of object x is used.
#' @param n.est.grid Number of points in the estimation grid, where hazard
#'   estimates are computed (used for \code{\link[muhaz]{muhaz}}). Default value
#'   is 101.
#' @param ask If \code{TRUE}, the user is asked for input, before a new figure
#'   is drawn.
#' @param ... Additional arguments passed through to plotting functions.
#' @seealso \code{\link[muhaz]{muhaz}}
#' @examples
#' data(survdata)
#' cp <- cpsurv(survdata$time, survdata$event, cpmax = 360, intwd = 10)
#' plot(cp, ask = FALSE)
#'
#' \dontrun{
#' cp <- cpsurv(survdata$time, survdata$event, cpmax = 360, intwd = 10,
#' boot.ci = TRUE)
#' plot(cp, type = "pvals", ask = FALSE)
#' }
#' @importFrom muhaz muhaz
#' @importFrom grDevices devAskNewPage
#' @importFrom graphics abline barplot lines plot
#' @method plot cpsurv
#' @export
plot.cpsurv <- function(x, type = "all", ci = TRUE, ci.type = c("perc", "norm"),
                        const.haz = TRUE, regline = TRUE, legend = TRUE,
                        xlim = NULL, ylim = NULL, main = NULL, xlab = NULL,
                        ylab = NULL, min.time, max.time,
                        n.est.grid = 101, ask = TRUE, ...){

  ci.type <- match.arg(ci.type)
  stopifnot(is.logical(ci), is.logical(const.haz), is.logical(regline))
  if(missing(min.time))
    min.time <- 0
  if(missing(max.time))
    max.time <- x$cpmax
  if(is.null(xlim)) xlim <- c(0, x$cpmax)


  cp <- ifelse(!is.null(x$cp.bc), x$cp.bc, x$cp)
  existci <- !is.null(x$ci.normal)

  if (ask) {
    oask <- devAskNewPage()
    devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  if (any(type == "all" | type == "pvals")){
    ylim1 <- c(0,1)
    if(is.null(xlab)){
      xlab1 <- "follow-up time"
    } else{xlab1 <- xlab}
    if(is.null(ylab)){
      ylab1 <- "p-values"
    } else{ylab1 <- ylab}

    plot(x$lower.lim, x$p.values, xlim = xlim, ylim = ylim1,
         xlab = xlab1, ylab = ylab1, main = main, ...)
    if(regline){
      lines(x = c(0, cp), y = c(0, 0), lwd = 2)
      lines(x = c(cp, x$lower.lim[length(x$lower.lim)]),
            y = c(x$pv.mean, x$pv.mean), lwd = 2)
      lines(x = c(cp, cp), y = c(0, x$pv.mean), lty = 3)
      if(legend)
        legend("topleft", "Regression line", col = 1, lwd = 2, bg = "white")
    }
  }


  if (any(type == "all" | type == "events")){
    if(is.null(ylab)){
      ylab1 <- "events per units at risk"
    } else{ylab1 <- ylab}
    if(is.null(main)){
      main1 <- "Relative frequency of events per interval"
    } else{main1 <- main}
    if(all(x$event == 1)){
      events <- table(x$event, cut(x$time, seq(xlim[1], xlim[2], by=x$intwd)))
    } else{
      events <- table(x$event, cut(x$time, seq(xlim[1], xlim[2], by=x$intwd)))[2,]
    }

    riskset <- (length(x$time) - c(0, cumsum(events)[-length(events)]))
    relevents <- events / riskset
    barplot(relevents, ylab = ylab1, xlab = xlab, main = main1,
            ylim = ylim, col = "gray")
  }


  if (any(type == "all" | type == "hazard")){
    mh <- muhaz::muhaz(times = x$time, delta = x$event, min.time = min.time,
                       max.time = max.time,n.est.grid = n.est.grid)
    if(is.null(xlab)){
      xlab1 <- "follow-up time"
    } else{xlab1 <- xlab}
    if(is.null(ylab)){
      ylab1 <- "hazard rate"
    } else{ylab1 <- ylab}
    plot(mh, xlim = xlim, main = main, xlab = xlab1, ylab = ylab1, ylim = ylim, ...)

    haztext <- NULL
    if(const.haz){
      exprate <- sum(x$event[x$time > x$cpmax] == 1) /
        sum(x$time[x$time > x$cpmax] - x$cpmax)
      abline(h = exprate, col = "grey60")
      haztext <- "const.hazard"
      lines(mh, mgp = c(2.5,1,0))
    }

    abline(v=cp)

    if(ci && existci){
      if(ci.type == "perc"){
        abline(v = x$ci.percent[1], lty = 2, col = "grey40")
        abline(v = x$ci.percent[2], lty = 2, col = "grey40")
        if(legend)
          legend("topright", c("changepoint", "perc.int", haztext),
                 col = c(1,"grey40","grey60"),lty = c(1,2,1), bg = "white")
      } else{
        abline(v = x$ci.normal[1], lty = 2, col = "grey40")
        abline(v = x$ci.normal[2], lty = 2, col = "grey40")
        if(legend)
          legend("topright", c("changepoint", "normal.ci", haztext),
                 col = c(1,"grey40","grey60"),lty = c(1,2,1), bg = "white")
      }
    } else if(legend){
      legend("topright", c("changepoint", haztext),
             col = c(1,"grey60"),lty = c(1,1), bg = "white")
    }
  }
}




