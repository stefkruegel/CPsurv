
test_that("estimation causes no errors", {
  data("survdata")
  rounded <- survdata
  rounded$time <- ceiling(rounded$time)
  type1cens <- survdata
  type1cens$event <- 1
  type1cens[type1cens$time >= 540, "event"] <- 0

  with(rounded, cpsurv(time = time, event = event, cpmax = 300, seed = 12345))

  with(survdata, cpsurv(time = time, event = event, cpmax = 300, biascorrect = FALSE,
                        censoring = "no"))

  with(survdata, cpsurv(time = time, event = event, cpmax = 300, cpmin = 15,
                        parametric = TRUE, seed = 12345))

  with(survdata, cpsurv(time = time, event = event, cpmax = 300, boot.ci = TRUE,
                        parallel = FALSE, B.correct = 10, B = 20, seed = 12345))

  with(rounded, cpsurv(time = time, event = event, cpmax = 300, boot.ci = TRUE,
                       parametric = TRUE, parallel = FALSE, B.correct = 10, B = 20,
                       seed = 12345))

  with(survdata, cpsurv(time = time, event = event, cpmax = 260, intwd = 20,
                        cpmin = 17, parametric = TRUE, seed = 12345, boot.ci = TRUE,
                        parallel = FALSE, B.correct = 10, B = 20))

  with(type1cens, cpsurv(time = time, event = event, cpmax = 300, intwd = 20,
                         seed = 12345, boot.ci = TRUE, censoring = "type1",
                         parallel = FALSE, B.correct = 10, B = 20, censpoint = 540))
})

