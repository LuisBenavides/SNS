#' Calibration of the control limit for the selected chart
#' @description The methodology used to calibrate the control limit
#' for the SNS chart depending on the selected chart
#' @inheritParams mgetARL
#' @param targetARL scalar. is the target ARL to calibrate. By default is set to NULL
#' @param targetMRL scalar. is the target ARL to calibrate. By default is set to NULL
#' @param maxIter scalar. is a numeric. The maximum number of iteration to take the calibration before stops
#' @note The argument \code{chart.par} in this function correspond to the initial parameters to start the calibration.
#' @export
mcalibrateControlLimit <- function(targetARL = NULL, targetMRL = NULL,
                                  n, m, nv, theta = NULL, Ftheta = NULL,
                                  dists=c("Normal", "Normal"), mu=c(0,0), sigma=NULL, dists.par = matrix(c(0,1,1,0,1,1), ncol=2), correlation=0,
                                  chart="T2", chart.par=c(10), replicates = 50000,
                                  isParallel = FALSE, maxIter = 20, progress = TRUE,
                                  alignment="unadjusted", constant=NULL, absolute=FALSE) {
  # Check for errors
  if (is.null(targetARL) && is.null(targetMRL)) {
    print("ERROR: Target ARL or target mRL missing")
    return()
  } else if (!is.null(targetARL) && !is.null(targetMRL)) {
    print("ERROR: Two targets defined, delete one")
    return()
  }
  p <- 0.1
  if (is.null(targetARL)) {
    ARL0 <- (targetMRL * 1.5) / 10
  } else {
    ARL0 <- targetARL
  }

  switch(chart,
         T2 = {
           name.par <- "h"
           index.par <- 1
         }
  )

  x <- rep(NA, maxIter)
  y <- x

  i <- 1
  x[i] <- chart.par[index.par]
  while (i < maxIter) {
    chart.par[index.par] <- x[i]
    result <- SNS::mgetARL(n = n, m = m, nv = nv,
                      theta = theta, Ftheta = Ftheta,
                      dists = dists, mu = mu, sigma = sigma, dists.par = dists.par,
                      correlation=correlation, chart = chart, chart.par = chart.par,
                      replicates = replicates, isParallel = isParallel, calibrate = TRUE, arl0 = targetARL,
                      alignment=alignment, constant=constant,absolute=absolute)
    if (!is.null(targetARL)) {
      y[i] <- result$ARL
      target <- targetARL
      name <- "ARL"
    } else {
      y[i] <- result$MRL
      target <- targetMRL
      name <- "MRL"
    }

    if (abs(y[i] - target) <= 0.05 * target) {
      if (progress) cat("Convergence found with", name.par, "=", x[i], "--", name, "=", y[i], "\n", sep = " ")
      output <- list(
        objective.function = y[i],
        par.value = x[i],
        found = TRUE
      )
      return(output)
    } else {
      f1 <- 0
      f2 <- 0
      if (i > 2) {
        f1 <- x[i] - target
        f2 <- x[i - 1] - target
      }

      if (f1 * f2 < 0) {
        x0 <- x[i - 1]
        x1 <- x[i]
        y0 <- y[i - 1]
        y1 <- y[i]
        m <- (y1 - y0) / (x1 - x0)
        b <- y0 - m * x0
        x2 <- (target - b) / m
        x[i + 1] <- x2
      } else {
        if (y[i] <= target) {
          x[i + 1] <- x[i] * (1 + p)
        } else {
          x[i + 1] <- x[i] * (1 - p)
        }
        if (progress) cat("obtained=", y[i], " target=", target, " Change h=", x[i], " to h=", x[i + 1], "\n", sep = "")
      }
    }
    i <- i + 1
  }

  posMin <- which.min(abs(target - y))
  if (progress) cat("Best", name.par, "found ", x[posMin], "--", name, "=", y[posMin], "\n", sep = " ")

  output <- list(
    objective.function = y[posMin],
    par.value = x[posMin],
    found = FALSE
  )
  return(output)
}
