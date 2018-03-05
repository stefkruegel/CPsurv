#library(testthat)

data("survdata")

test_that("input checking works", {
  datatrunc <- survdata
  datatrunc[datatrunc$time > 360, "event"] <- 0
  expect_error(with(datatrunc, cpsurv(time, event, intwd = 20, cpmax = 360)),
               "No events with 'time' > 'cpmax'")

  wrongevent <- survdata$event
  wrongevent[5] <- 0.5
  expect_error(cpsurv(survdata$time, wrongevent, intwd = 20, cpmax = 360),
               "Argument 'event' has to be binary")

  expect_error(cpsurv(survdata$time, survdata$event[-1], intwd = 20,
                      cpmax = 360),
               "Vectors 'time' and 'event' must be of equal length.")

  expect_error(with(survdata, cpsurv(time, event, intwd = "nonsense",
                                     cpmax = 360)),
               "Argument 'intwd' is not a single numeric value")

  expect_error(with(survdata, cpsurv(time, event, conf.level = 2,
                                     cpmax = 360)),
               "Value for argument 'conf.level' too high")
})
