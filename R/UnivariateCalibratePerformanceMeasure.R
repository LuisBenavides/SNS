#' Calibration of the control limit for the selected chart
#' @description The methodology used to calibrate the control limit
#' for the SNS chart depending on the selected chart
#' @inheritParams getARL
#' @param targetARL scalar. is the target ARL to calibrate. By default is set to NULL
#' @param targetMRL scalar. is the target ARL to calibrate. By default is set to NULL
#' @param maxIter scalar. is a numeric. The maximum number of iteration to take the calibration before stops
#' @note The argument \code{chart.par} in this function correspond to the initial parameters to start the calibration.
#' @export
#' @examples
#' n <- 2 # subgroup size
#' m <- 30 # reference-sample size
#' dist <- "Normal" # distribution
#' mu <- c(0, 0) # c(reference sample mean, monitoring sample mean)
#' sigma <- c(1, 1) # c(reference sample sd, monitoring sample sd)
#'
#' #### Distribution parameters
#' dist.par <- c(0, 1) # c(location, scale)
#'
#' #### Other Parameters
#' replicates <- 2
#' targetARL <- 370
#' isParallel = FALSE
#'
#' #### Control chart parameters
#' chart <- "Shewhart"
#' chart.par <- c(3)
#' shewhart <- calibrateControlLimit(
#'   targetARL = targetARL, targetMRL = NULL, n = n, m = m, theta = NULL,
#'   Ftheta = NULL, dist = dist, mu = mu, sigma = sigma, dist.par = dist.par, chart.par = chart.par,
#'   replicates = replicates, chart = chart, isParallel = isParallel
#' )
#'
#' chart <- "CUSUM"
#' chart.par <- c(0.5, 2.5, 3)
#' cusum <- calibrateControlLimit(
#'   targetARL = targetARL, targetMRL = NULL, n = n, m = m, theta = NULL,
#'   Ftheta = NULL, dist = dist, mu = mu, sigma = sigma, dist.par = dist.par, chart.par = chart.par,
#'   replicates = replicates, chart = chart, isParallel = isParallel
#' )
#'
#' chart <- "EWMA"
#' chart.par <- c(0.2, 2.962)
#' ewma <- calibrateControlLimit(
#'   targetARL = targetARL, targetMRL = NULL, n = n, m = m, theta = NULL,
#'   Ftheta = NULL, dist = dist, mu = mu, sigma = sigma, dist.par = dist.par, chart.par = chart.par,
#'   replicates = replicates, chart = chart, isParallel = isParallel
#' )
calibrateControlLimit <- function(targetARL = NULL, targetMRL = NULL,
                                  n, m, theta = NULL, Ftheta = NULL,
                                  scoring = "Z", Chi2corrector="None",
                                  dist, mu, sigma, dist.par = c(0, 1, 1),
                                  chart, chart.par, replicates = 50000,
                                  isParallel = TRUE, maxIter = 20, progress = TRUE,
                                  alignment="unadjusted", constant=NULL, absolute=FALSE,
                                  isFixed=FALSE, rounding.factor = NULL) {
  # Check for errors
  if (is.null(targetARL) && is.null(targetMRL)) {
    print("ERROR: Target ARL or target mRL missing")
    return()
  } else if (!is.null(targetARL) && !is.null(targetMRL)) {
    print("ERROR: Two targets defined, delete one")
    return()
  }
  #auxiliar variable to control when the interpolation start
  startInterpolate = FALSE
  #percentage of increment for par.value
  p <- 0.1
  if (is.null(targetARL)) {
    # if MRL is selected
    # set the maximum ARL in getARL function to
    # at least the expected ARL for that MRL around MRL/0.7
    ARL0 <- targetMRL / 0.7
  }else {
    # if MRL is selected
    # set the maximum ARL in getARL function to
    # at 100 times the ARL
    ARL0 <- targetARL * 100
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

  x <- rep(NA, maxIter)
  y <- x

  i <- 1
  x[i] <- chart.par[index.par]


  while (i < maxIter) {
    chart.par[index.par] <- x[i]

    result <- SNS::getARL(n = n, m = m, theta = theta, Ftheta = Ftheta,
                          dist = dist, mu = mu, sigma = sigma, dist.par = dist.par,
                          chart = chart, chart.par = chart.par, replicates = replicates,
                          isParallel = isParallel, calibrate = TRUE, arl0 = ARL0,
                          alignment=alignment, constant=constant,absolute=absolute,isFixed=isFixed,
                          scoring=scoring,Chi2corrector=Chi2corrector, rounding.factor = rounding.factor)
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
      #if the obtained value is in its 5% from target
      #return the par.value and the obtained value
      if (progress) cat("Convergence found with", name.par, "=", x[i], "--", name, "=", y[i], "\n", sep = " ")
      output <- list(
        objective.function = y[i],
        par.value = x[i],
        iter = i,
        found = TRUE
      )
      return(output)
    } else {
      fi <- y[i] - target
      if (!startInterpolate){
        #if only increase or decrease values (not enter to interpolation)
        #update values point 0
        f0 <- y[i-1] - target
        x0 <- x[i-1]
        y0 <- y[i-1]

        #update values point 1
        f1 <- y[i] - target
        x1 <- x[i]
        y1 <- y[i]
      }
      if (i >= 2){#when at least are two values
        #obtain error
        f0 <- y0 - target
        f1 <- y1 - target
      }else{#when is the first iteration (initilize values)
        f0 <- 0
        f1 <- 0
      }
      if(fi * f0 < 0 && fi != f0){
        #if there is a sign change considering f0 as the comparison
        #fi != f0 needed for first iteration and repeated vaules
        startInterpolate = TRUE
        #update point (x1,y1,f1)
        x1 <- x[i]
        y1 <- y[i]
        f1 <- fi

        #interpolate
        m1 <- (y1 - y0) / (x1 - x0)
        b <- y0 - m1 * x0
        x2 <- (target - b) / m1
        x[i + 1] <- x2
      }else if(fi * f1 < 0 && fi != f1){
        #if there is a sign change considering f1 as the comparison
        #fi != f0 needed for first iteration and repeated vaules
        startInterpolate = TRUE

        #update point (x0,y0,f0)
        x0 <- x[i]
        y0 <- y[i]
        f0 <- fi

        #interpolate
        m1 <- (y1 - y0) / (x1 - x0)
        b <- y0 - m1 * x0
        x2 <- (target - b) / m1
        x[i + 1] <- x2
      }else {
        if (y[i] <= target) {
          #if the target is not reached and the obtained value is below
          #increase p percent its value
          x[i + 1] <- x[i] * (1 + p)
        } else {
          #if the target is not reached and the obtained value is above
          #decrease p percent its value
          x[i + 1] <- x[i] * (1 - p)
        }
        if (progress) cat("obtained=", y[i], " target=", target, " Change h=", x[i], " to h=", x[i + 1], "\n", sep = "")
      }
    }
    i <- i + 1
  }

  #if run out of itertations (target value not encountered)
  #get the par.value that has the closest value to the target
  posMin <- which.min(abs(target - y))
  if (progress) cat("Best", name.par, "found ", x[posMin], "--", name, "=", y[posMin], "\n", sep = " ")

  output <- list(
    objective.function = y[posMin],
    par.value = x[posMin],
    iter = i,
    found = FALSE
  )
  return(output)
}

