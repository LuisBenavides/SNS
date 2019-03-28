#' Calibration of the control limit for the selected chart
#'
#' The methodology used to calibrate the control limit
#' for the SNS chart depending on the chart selected
#'
#' @param targetARL is the target ARL to calibrate. By default is set to NULL
#' @param targetMRL is the target ARL to calibrate. By default is set to NULL
#' @param n is the subroup size
#' @param m is the reference sample size
#' @param theta is a constant. By default is set to NULL
#' @param Ftheta is a constant between (0,1). By default is set to NULL
#' @param dist is the distribution to analyze. By default is "Normal"
#' @param mu is a vector of two elements. By default is set to c(0,0)
#' The firt element refers to the mean of the reference sample.
#' The second element refers to the mean of the monitoring sample.
#' @param sigma is a vector of two elements. By default is set to c(1,1)
#' The firt element refers to the standard deviation of the reference sample.
#' The second element refers to the standard deviation of the monitoring sample.
#' @param dist.par is a vector of three elements. By defautl is set to c(0, 1, 1)
#' The firt element refers to location of the distribution
#' The second element refers to scale of the distribution
#' The third element refers to shape of the distribution
#' @param chart is the selected type of chart. Three options are available: Shewhart, CUSUM, EWMA
#' @param initial.par is a vector of one, or two elements depending
#' on the chart selected: Shewhart, CUSUM, EWMA
#' Shewhart chart is c(k), where k comes from UCL = mu + k*sigma, LCL = mu - k*sigma.
#' CUSUM chart is c(k, h, t) where k is the reference value and h is the control limit,
#' and t is the type of the chart (1:positive, 2:negative, 3:two sides)
#' EWMA chart is c(lambda, L), where lambda is the smoothing constant
#' and L multiplies standard deviation to get the control limit
#' @param replicates is a numeric and represent the number of replicates
#' to run the function getARL in each iteration
#' @param isParallel is a boolean. If is TRUE the code runs in parallel according to the
#' number of cores in the computer,otherwise the code runs sequentially. By default is set to TRUE
#' @param maxIter is a numeric. The maximum number of iteration to take the calibration before stops
#' @param progress is a boolean. If TRUE it shows the progress of the calibration in the console.
#' @export
#' @examples
#' n <- 2 # subgroup size
#' m <- 30 # reference-sample size
#' dist <- "Normal" # distribution
#' mu <- c(0, 0) # c(reference sample mean, monitoring sample mean)
#' sigma <- c(1, 1) # c(reference sample sd, monitoring sample sd)
#'
#' #### Distribution parameters
#' dist.par <- c(0, 1, 1) # c(location, scale, shape)
#'
#' #### Other Parameters
#' replicates <- 2
#' targetARL <- 370
#'
#' #### Control chart parameters
#' chart <- "Shewhart"
#' chart.par <- c(3)
#' shewhart <- calibrateControlLimit(
#'   targetARL = targetARL, targetMRL = NULL, n = n, m = m, theta = NULL,
#'   Ftheta = NULL, dist = dist, mu = mu, sigma = sigma, dist.par = dist.par, initial.par = chart.par,
#'   replicates = replicates, chart = chart
#' )
#'
#' chart <- "CUSUM"
#' chart.par <- c(0.5, 2.5, 3)
#' cusum <- calibrateControlLimit(
#'   targetARL = targetARL, targetMRL = NULL, n = n, m = m, theta = NULL,
#'   Ftheta = NULL, dist = dist, mu = mu, sigma = sigma, dist.par = dist.par, initial.par = chart.par,
#'   replicates = replicates, chart = chart
#' )
#'
#' chart <- "EWMA"
#' chart.par <- c(0.2, 2.962)
#' ewma <- calibrateControlLimit(
#'   targetARL = targetARL, targetMRL = NULL, n = n, m = m, theta = NULL,
#'   Ftheta = NULL, dist = dist, mu = mu, sigma = sigma, dist.par = dist.par, initial.par = chart.par,
#'   replicates = replicates, chart = chart
#' )
calibrateControlLimit <- function(targetARL = NULL, targetMRL = NULL, n, m, theta = NULL, Ftheta = NULL, dist, mu, sigma, dist.par = c(0, 1, 1), chart, initial.par, replicates = 50000, isParallel = FALSE, maxIter = 20, progress = TRUE) {
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
    Shewhart = {
      name.par <- "k"
      index.par <- 1
    },
    CUSUM = {
      name.par <- "h"
      index.par <- 2
    },
    EWMA = {
      name.par <- "L"
      index.par <- 2
    }
  )

  x <- rep(NA, 3)
  y <- x

  i <- 1
  x[i] <- initial.par[index.par]
  while (i < maxIter) {
    chart.par[index.par] <- x[i]
    result <- getARL(n = n, m = m, theta = theta, Ftheta = Ftheta, dist = dist, mu = mu, sigma = sigma, dist.par = dist.par, chart = chart, chart.par = chart.par, replicates = replicates, isParallel = isParallel, calibrate = TRUE, arl0 = targetARL)
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
