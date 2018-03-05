context("test methods")

data(survdata)
cpest1 <- cpsurv(time = survdata$time, event = survdata$event, cpmax = 300,
                 boot.ci = FALSE, seed = 12345)
cpest2 <- cpsurv(time = survdata$time, event = survdata$event, cpmax = 300,
                 biascorrect = TRUE, boot.ci = FALSE, seed = 12345)
cpest3 <- cpsurv(time = survdata$time, event = survdata$event, cpmax = 300,
                 biascorrect = TRUE, boot.ci = TRUE, parallel = FALSE, B = 5, seed = 12345)
cpest4 <- cpsurv(time = survdata$time, event = survdata$event, cpmax = 300,
                 cpmin = 20, intwd = 10, boot.ci = TRUE, parallel = FALSE, B = 5, seed = 12345)

test_that("print doesn't cause errors", {
  print(cpest1)
  print(cpest2)
  print(cpest3)
  print(cpest4)
})

test_that("summary doesn't cause errors", {
  summary(cpest1)
  summary(cpest2)
  summary(cpest3)
  summary(cpest4)
})

test_that("plot has right input checking", {
  expect_error(plot(cpest3, ci = "nonsense"))
})

test_that("plot works", {
  plot(cpest3)
  plot(cpest4, legend = FALSE)
  plot(cpest3, ci.type = "norm")
  plot(cpest3, type = "hazard")
  plot(cpest3, type = c("pval", "events"))
  plot(cpest3, ci = FALSE, const.haz = FALSE, regline = F, xlim = c(50, 200))
})

